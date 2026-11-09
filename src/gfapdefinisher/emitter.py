"""Sequence and AGP assembly emitter with empirical provenance tracking.

Stitches unitig paths, rDNA cassettes, and bridged gaps into finished FASTA
and AGP v2.1 outputs with empirical provenance metadata and unpolished sequence flags.
Implements BIO-004 (repeat preservation) and BIO-005 (lossless overlap trimming).
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any, TextIO

from gfapdefinisher.models import Strand
from gfapdefinisher.policy import DEFAULT_POLICY, Policy
from gfapdefinisher.policy.models import PolicyValue
from gfapdefinisher.telomere import scan_terminal_tract

__all__ = [
    "AgpGapRow",
    "AgpPolicy",
    "AgpRecord",
    "BreakSiteGenomeSummary",
    "DEFAULT_AGP_POLICY",
    "FastaAgpEmitter",
    "PathComponent",
    "ScaffoldPlan",
    "ScaffoldStep",
    "SequenceProvenance",
    "SequenceSummaryMetrics",
    "assemble_path_sequence",
    "assemble_scaffold_sequence",
    "compute_break_site_summary",
    "emit_agp",
    "emit_excluded_fasta",
    "emit_fasta",
    "emit_grs_bed",
    "emit_rdna_bed",
    "emit_scaffold",
    "emit_summary_tsv",
    "emit_telomere_bed",
    "find_break_site_scaffolds",
    "format_summary_table",
    "nominate_super_scaffolds",
    "parse_agp",
    "parse_summary_tsv",
]


@dataclass(frozen=True)
class AgpPolicy:
    """AGP emission policy parameters."""

    unknown_gap_bp: PolicyValue[int] = PolicyValue(
        value=int("100"),
        justification="Standard AGP v2.1 gap length for unknown gaps (type U, repeat/tangle).",
        provenance="NCBI AGP specification v2.1 standard convention",
    )


DEFAULT_AGP_POLICY = AgpPolicy()


class SequenceProvenance(str, Enum):
    """Empirical sequence provenance and risk categorization for assembled components."""

    HIFI_CORE = "HiFi-core"
    HIFI_GERMLINE = "HiFi-germline"
    HIFI_1X = "HiFi-1x"
    ONT_DERIVED = "ONT-derived"
    RDNA_CHAINED = "rDNA-chained"

    # Lowercase aliases
    hifi_core = "HiFi-core"
    hifi_germline = "HiFi-germline"
    hifi_1x = "HiFi-1x"
    ont_derived = "ONT-derived"
    rdna_chained = "rDNA-chained"


@dataclass(slots=True)
class PathComponent:
    """Component of an assembled scaffold path.

    Attributes:
        segment_id: Assembly graph unitig or read identifier.
        strand: Orientation along scaffold ('+' or '-').
        sequence: Oriented nucleotide sequence of this component.
        overlap_to_next_bp: Exact dovetail overlap length to next component in bp.
        sequence_provenance: Empirical provenance category.
        unpolished: Whether sequence is unpolished raw long-read data.
        is_gap: Whether this component represents an unsequenced AGP gap.
        gap_length: Length of gap in bp.
        gap_type: AGP gap type (e.g. 'repeat', 'scaffold').
        linkage: AGP linkage flag ('yes' or 'no').
        linkage_evidence: AGP linkage evidence (e.g. 'rna_tangle', 'hifi-spanning').
        gap_component_type: AGP component designation ('N' or 'U').
    """

    segment_id: str = ""
    strand: Strand | str = "+"
    sequence: str = ""
    overlap_to_next_bp: int = 0
    sequence_provenance: SequenceProvenance = SequenceProvenance.HIFI_CORE
    unpolished: bool = False
    is_gap: bool = False
    gap_length: int = 0
    gap_type: str = "scaffold"
    linkage: str = "yes"
    linkage_evidence: str = "hifi-spanning"
    gap_component_type: str = "N"
    feature_id: str | None = None

    def __init__(
        self,
        segment_id: str | None = None,
        strand: Strand | str = "+",
        sequence: str = "",
        overlap_to_next_bp: int = 0,
        sequence_provenance: SequenceProvenance | str = SequenceProvenance.HIFI_CORE,
        unpolished: bool = False,
        *,
        component_id: str | None = None,
        overlap_to_next: int | None = None,
        provenance: SequenceProvenance | str | None = None,
        is_gap: bool = False,
        gap_length: int = 0,
        gap_type: str = "scaffold",
        linkage: str = "yes",
        linkage_evidence: str = "hifi-spanning",
        gap_component_type: str = "N",
        feature_id: str | None = None,
    ) -> None:
        seg_id = segment_id if segment_id is not None else (component_id or "")
        ov = overlap_to_next if overlap_to_next is not None else overlap_to_next_bp
        prov = provenance if provenance is not None else sequence_provenance
        if is_gap and not sequence and gap_length > 0:
            sequence = "N" * gap_length

        if not seg_id and not is_gap:
            raise ValueError("segment_id cannot be empty.")
        if not seg_id and is_gap:
            seg_id = "gap"

        strand_val = strand.value if isinstance(strand, Strand) else strand
        if strand_val not in ("+", "-"):
            raise ValueError(f"Invalid strand: {strand!r}. Must be '+' or '-'.")

        if ov < 0:
            raise ValueError(f"overlap_to_next_bp cannot be negative: {ov}")
        if not is_gap and ov > len(sequence):
            raise ValueError(f"overlap_to_next_bp ({ov}) cannot exceed sequence length ({len(sequence)}).")

        if isinstance(prov, str) and not isinstance(prov, SequenceProvenance):
            prov = SequenceProvenance(prov)

        object.__setattr__(self, "segment_id", seg_id)
        object.__setattr__(self, "strand", strand_val)
        object.__setattr__(self, "sequence", sequence)
        object.__setattr__(self, "overlap_to_next_bp", ov)
        object.__setattr__(self, "sequence_provenance", prov)
        object.__setattr__(self, "unpolished", unpolished)
        object.__setattr__(self, "is_gap", is_gap)
        object.__setattr__(self, "gap_length", gap_length if gap_length > 0 else len(sequence))
        object.__setattr__(self, "gap_type", gap_type)
        object.__setattr__(self, "linkage", linkage)
        object.__setattr__(self, "linkage_evidence", linkage_evidence)
        object.__setattr__(self, "gap_component_type", gap_component_type)
        object.__setattr__(self, "feature_id", feature_id)

    @property
    def is_unpolished(self) -> bool:
        """Alias for unpolished flag."""
        return self.unpolished


@dataclass(slots=True)
class AgpRecord:
    """AGP v2.1 component record with empirical provenance metadata.

    Attributes:
        object_id: Scaffold or chromosome object identifier.
        object_beg: 1-based start coordinate in the object.
        object_end: 1-based end coordinate in the object.
        part_number: 1-based order of this component in the object.
        component_type: AGP component type (e.g. 'W' for sequence, 'N' or 'U' for gap).
        component_id: Component identifier or gap length.
        component_beg: 1-based start coordinate on component, or gap type.
        component_end: 1-based end coordinate on component, or linkage.
        orientation: Orientation ('+', '-') or linkage evidence.
        sequence_provenance: Empirical provenance category.
        unpolished: Whether component sequence is unpolished.
        comment: Optional auxiliary comment text.
    """

    object_id: str
    object_beg: int
    object_end: int
    part_number: int
    component_type: str
    component_id: str
    component_beg: int | str
    component_end: int | str
    orientation: str
    sequence_provenance: SequenceProvenance | str | None = None
    unpolished: bool = False
    comment: str | None = None

    def __post_init__(self) -> None:
        if not self.object_id:
            raise ValueError("object_id cannot be empty.")
        if self.object_beg < 1:
            raise ValueError(f"object_beg must be >= 1, got {self.object_beg}")
        if self.object_end < self.object_beg:
            raise ValueError(f"object_end ({self.object_end}) cannot be less than object_beg ({self.object_beg})")
        if self.part_number < 1:
            raise ValueError(f"part_number must be >= 1, got {self.part_number}")
        if isinstance(self.sequence_provenance, str) and not isinstance(self.sequence_provenance, SequenceProvenance):
            try:
                object.__setattr__(self, "sequence_provenance", SequenceProvenance(self.sequence_provenance))
            except ValueError:
                pass

    @property
    def is_unpolished(self) -> bool:
        """Alias for unpolished flag."""
        return self.unpolished

    def to_tsv_row(self) -> list[str]:
        """Serialize AGP record to ordered list of column strings."""
        row = [
            self.object_id,
            str(self.object_beg),
            str(self.object_end),
            str(self.part_number),
            self.component_type,
            str(self.component_id),
            str(self.component_beg),
            str(self.component_end),
            self.orientation,
        ]
        if self.sequence_provenance is not None:
            prov_str = (
                self.sequence_provenance.value
                if isinstance(self.sequence_provenance, Enum)
                else str(self.sequence_provenance)
            )
            tag = f"# provenance={prov_str};unpolished={self.unpolished}"
            if self.comment:
                tag = f"{tag}; {self.comment}"
            row.append(tag)
        elif self.comment:
            row.append(f"# {self.comment}")
        return row

    def to_line(self) -> str:
        """Format as a tab-delimited AGP v2.1 row with optional provenance comment."""
        return "\t".join(self.to_tsv_row())


def assemble_scaffold_sequence(
    components: Sequence[PathComponent],
    scaffold_id: str = "scaffold",
) -> tuple[str, list[AgpRecord]]:
    """Assemble oriented path components into a stitched scaffold sequence and AGP records.

    Stitches oriented component sequences while trimming exact dovetail overlaps
    between adjacent components according to overlap_to_next_bp (BIO-005).

    Args:
        components: Sequence of PathComponent objects along the scaffold path.
        scaffold_id: Identifier for the assembled scaffold/chromosome object.

    Returns:
        Tuple of (assembled_nucleotide_sequence, list_of_agp_records).
    """
    if not components:
        return "", []

    current_scaffold_pos = 1
    part_number = 1
    chunks: list[str] = []
    records: list[AgpRecord | AgpGapRow] = []
    prev_overlap = 0

    for comp in components:
        if comp.is_gap:
            gap_len = comp.gap_length if comp.gap_length > 0 else len(comp.sequence)
            chunks.append("N" * gap_len)
            obj_beg = current_scaffold_pos
            obj_end = current_scaffold_pos + gap_len - 1
            current_scaffold_pos = obj_end + 1
            comp_type = comp.gap_component_type
            if (
                comp_type == "N"
                and gap_len == DEFAULT_AGP_POLICY.unknown_gap_bp.value
                and comp.linkage_evidence in ("rna_tangle", "hifi-spanning")
            ):
                comp_type = "U"
            row = AgpGapRow(
                object_id=scaffold_id,
                object_beg=obj_beg,
                object_end=obj_end,
                part_number=part_number,
                component_type=comp_type,
                gap_length=gap_len,
                gap_type=comp.gap_type,
                linkage=comp.linkage,
                linkage_evidence=comp.linkage_evidence,
            )
            records.append(row)
            part_number += 1
            prev_overlap = 0
            continue

        strand_str = comp.strand if isinstance(comp.strand, str) else comp.strand.value
        seq_len = len(comp.sequence)
        prefix_trim = prev_overlap

        if prefix_trim >= seq_len:
            raise ValueError(
                f"Overlap {prefix_trim} bp exceeds or equals component sequence length {seq_len} bp for {comp.segment_id}"
            )

        trimmed_seq = comp.sequence[prefix_trim:]
        chunk_len = len(trimmed_seq)
        chunks.append(trimmed_seq)

        obj_beg = current_scaffold_pos
        obj_end = current_scaffold_pos + chunk_len - 1
        current_scaffold_pos = obj_end + 1

        if strand_str == "+":
            comp_beg = 1 + prefix_trim
            comp_end = seq_len
        else:
            comp_beg = 1
            comp_end = seq_len - prefix_trim

        rec = AgpRecord(
            object_id=scaffold_id,
            object_beg=obj_beg,
            object_end=obj_end,
            part_number=part_number,
            component_type="W",
            component_id=comp.segment_id,
            component_beg=comp_beg,
            component_end=comp_end,
            orientation=strand_str,
            sequence_provenance=comp.sequence_provenance,
            unpolished=comp.unpolished,
        )
        records.append(rec)
        part_number += 1
        prev_overlap = comp.overlap_to_next_bp

    assembled_seq = "".join(chunks)
    return assembled_seq, records


# Semantic alias
assemble_path_sequence = assemble_scaffold_sequence


@dataclass(slots=True, init=False)
class ScaffoldStep:
    """Step in an assembled scaffold plan representing a segment, ONT gap, or unknown gap."""

    segment_id: str
    strand: Strand | str
    sequence: str | None
    overlap_to_next_bp: int
    provenance: SequenceProvenance
    unpolished: bool
    gap_resolution: Any
    is_unknown_gap: bool
    unknown_gap_type: str
    unknown_gap_evidence: str

    def __init__(
        self,
        segment_id: str = "",
        strand: Strand | str = Strand.FORWARD,
        sequence: str | None = None,
        overlap_to_next_bp: int = 0,
        provenance: SequenceProvenance | str | None = None,
        sequence_provenance: SequenceProvenance | str | None = None,
        unpolished: bool = False,
        gap_resolution: Any = None,
        is_unknown_gap: bool = False,
        unknown_gap_type: str = "repeat",
        unknown_gap_evidence: str = "rna_tangle",
    ) -> None:
        if isinstance(strand, Strand):
            strand_val = strand.value
        elif strand in ("+", "-"):
            strand_val = strand
        else:
            raise ValueError(f"Invalid strand: {strand!r}. Must be '+' or '-'.")

        prov = provenance if provenance is not None else sequence_provenance
        if prov is None:
            prov_val = SequenceProvenance.HIFI_CORE
        elif isinstance(prov, str) and not isinstance(prov, SequenceProvenance):
            prov_val = SequenceProvenance(prov)
        else:
            prov_val = prov

        object.__setattr__(self, "segment_id", segment_id)
        object.__setattr__(self, "strand", strand_val)
        object.__setattr__(self, "sequence", sequence)
        object.__setattr__(self, "overlap_to_next_bp", overlap_to_next_bp)
        object.__setattr__(self, "provenance", prov_val)
        object.__setattr__(self, "unpolished", unpolished)
        object.__setattr__(self, "gap_resolution", gap_resolution)
        object.__setattr__(self, "is_unknown_gap", is_unknown_gap)
        object.__setattr__(self, "unknown_gap_type", unknown_gap_type)
        object.__setattr__(self, "unknown_gap_evidence", unknown_gap_evidence)

    @property
    def sequence_provenance(self) -> SequenceProvenance:
        """Alias for provenance."""
        return self.provenance

    @property
    def is_unpolished(self) -> bool:
        """Alias for unpolished flag."""
        return self.unpolished


@dataclass(slots=True)
class ScaffoldPlan:
    """Plan for assembling a scaffold from ordered ScaffoldSteps."""

    scaffold_id: str
    steps: list[ScaffoldStep] = field(default_factory=list)


@dataclass(slots=True)
class AgpGapRow:
    """AGP v2.1 gap component record."""

    object_id: str
    object_beg: int
    object_end: int
    part_number: int
    component_type: str
    gap_length: int
    gap_type: str
    linkage: str
    linkage_evidence: str
    comment: str | None = None

    def __post_init__(self) -> None:
        if not self.object_id:
            raise ValueError("object_id cannot be empty.")
        if self.object_beg < 1:
            raise ValueError(f"object_beg must be >= 1, got {self.object_beg}")
        if self.object_end < self.object_beg:
            raise ValueError(f"object_end ({self.object_end}) cannot be less than object_beg ({self.object_beg})")
        if self.part_number < 1:
            raise ValueError(f"part_number must be >= 1, got {self.part_number}")

    @property
    def component_id(self) -> str:
        return str(self.gap_length)

    @property
    def component_beg(self) -> str:
        return self.gap_type

    @property
    def component_end(self) -> str:
        return self.linkage

    @property
    def orientation(self) -> str:
        return self.linkage_evidence

    @property
    def sequence_provenance(self) -> None:
        return None

    @property
    def unpolished(self) -> bool:
        return False

    @property
    def is_unpolished(self) -> bool:
        return False

    def to_tsv_row(self) -> list[str]:
        row = [
            self.object_id,
            str(self.object_beg),
            str(self.object_end),
            str(self.part_number),
            self.component_type,
            str(self.gap_length),
            self.gap_type,
            self.linkage,
            self.linkage_evidence,
        ]
        if self.comment:
            row.append(f"# {self.comment}")
        return row

    def to_line(self) -> str:
        return "\t".join(self.to_tsv_row())


class FastaAgpEmitter:
    """FASTA and AGP v2.1 scaffold sequence and record emitter."""

    def __init__(
        self,
        policy: Policy = DEFAULT_POLICY,
        sequences: Mapping[str, str] | None = None,
    ) -> None:
        self.policy = policy
        if not hasattr(self.policy, "agp"):
            object.__setattr__(self.policy, "agp", AgpPolicy())
        self.sequences = sequences or {}

    def emit_scaffold(
        self,
        plan: ScaffoldPlan,
        fasta_output: Path | TextIO | str | None = None,
        agp_output: Path | TextIO | str | None = None,
    ) -> tuple[str, list[AgpRecord | AgpGapRow]]:
        """Assemble a ScaffoldPlan into contiguous sequence and AGP records."""
        current_pos = 1
        part_num = 1
        seq_parts: list[str] = []
        rows: list[AgpRecord | AgpGapRow] = []
        prev_overlap = 0

        for step in plan.steps:
            if step.is_unknown_gap:
                gap_len = self.policy.agp.unknown_gap_bp.value
                seq_parts.append("N" * gap_len)
                row = AgpGapRow(
                    object_id=plan.scaffold_id,
                    object_beg=current_pos,
                    object_end=current_pos + gap_len - 1,
                    part_number=part_num,
                    component_type="U",
                    gap_length=gap_len,
                    gap_type=step.unknown_gap_type,
                    linkage="yes",
                    linkage_evidence=step.unknown_gap_evidence,
                )
                rows.append(row)
                current_pos += gap_len
                part_num += 1
                prev_overlap = 0
            elif step.gap_resolution is not None:
                gap_res = step.gap_resolution
                rep_seq = getattr(gap_res, "representative_sequence", None)
                if rep_seq:
                    chunk_len = len(rep_seq)
                    seq_parts.append(rep_seq)
                    obj_beg = current_pos
                    obj_end = current_pos + chunk_len - 1
                    current_pos = obj_end + 1
                    rep_read = getattr(gap_res, "representative_read", None)
                    read_id = getattr(rep_read, "qname", "ont_bridge") if rep_read else "ont_bridge"
                    rec = AgpRecord(
                        object_id=plan.scaffold_id,
                        object_beg=obj_beg,
                        object_end=obj_end,
                        part_number=part_num,
                        component_type="W",
                        component_id=read_id,
                        component_beg=1,
                        component_end=chunk_len,
                        orientation="+",
                        sequence_provenance=SequenceProvenance.ONT_DERIVED,
                        unpolished=True,
                    )
                    rows.append(rec)
                    part_num += 1
                    prev_overlap = 0
                else:
                    gap_len = self.policy.agp.unknown_gap_bp.value
                    seq_parts.append("N" * gap_len)
                    row = AgpGapRow(
                        object_id=plan.scaffold_id,
                        object_beg=current_pos,
                        object_end=current_pos + gap_len - 1,
                        part_number=part_num,
                        component_type="U",
                        gap_length=gap_len,
                        gap_type="scaffold",
                        linkage="yes",
                        linkage_evidence="unbridged",
                    )
                    rows.append(row)
                    current_pos += gap_len
                    part_num += 1
                    prev_overlap = 0
            elif step.segment_id:
                raw_seq = step.sequence
                if raw_seq is None and step.segment_id in self.sequences:
                    raw_seq = self.sequences[step.segment_id]
                if raw_seq is None:
                    raw_seq = ""

                strand_str = step.strand if isinstance(step.strand, str) else step.strand.value
                seq_len = len(raw_seq)
                prefix_trim = prev_overlap

                if prefix_trim > 0 and prefix_trim >= seq_len:
                    raise ValueError(
                        f"Overlap {prefix_trim} bp exceeds or equals component sequence length {seq_len} bp for {step.segment_id}"
                    )

                trimmed_seq = raw_seq[prefix_trim:]
                chunk_len = len(trimmed_seq)
                seq_parts.append(trimmed_seq)

                obj_beg = current_pos
                obj_end = current_pos + chunk_len - 1
                current_pos = obj_end + 1

                if strand_str == "+":
                    comp_beg = 1 + prefix_trim
                    comp_end = seq_len
                else:
                    comp_beg = 1
                    comp_end = seq_len - prefix_trim

                rec = AgpRecord(
                    object_id=plan.scaffold_id,
                    object_beg=obj_beg,
                    object_end=obj_end,
                    part_number=part_num,
                    component_type="W",
                    component_id=step.segment_id,
                    component_beg=comp_beg,
                    component_end=comp_end,
                    orientation=strand_str,
                    sequence_provenance=step.provenance,
                    unpolished=step.unpolished,
                )
                rows.append(rec)
                part_num += 1
                prev_overlap = step.overlap_to_next_bp

        assembled_seq = "".join(seq_parts)
        if fasta_output is not None:
            emit_fasta({plan.scaffold_id: assembled_seq}, fasta_output)
        if agp_output is not None:
            emit_agp(rows, agp_output)

        return assembled_seq, rows


def emit_scaffold(
    plan: ScaffoldPlan,
    policy: Policy = DEFAULT_POLICY,
    fasta_output: Path | TextIO | str | None = None,
    agp_output: Path | TextIO | str | None = None,
) -> tuple[str, list[AgpRecord | AgpGapRow]]:
    """Convenience function to emit a ScaffoldPlan using FastaAgpEmitter."""
    return FastaAgpEmitter(policy=policy).emit_scaffold(
        plan=plan,
        fasta_output=fasta_output,
        agp_output=agp_output,
    )


def emit_fasta(
    scaffolds: Mapping[str, str],
    output: Path | TextIO | str,
    line_width: int | None = None,
) -> None:
    """Emit scaffold sequences to a FASTA format destination.

    Args:
        scaffolds: Mapping from scaffold/chromosome identifier to sequence.
        output: Destination file path or open text stream.
        line_width: Optional line wrap width (None for unwrapped).
    """

    def _write(writer: TextIO) -> None:
        for scaffold_id, seq in scaffolds.items():
            writer.write(f">{scaffold_id}\n")
            if line_width is not None and line_width > 0:
                for pos in range(0, len(seq), line_width):
                    writer.write(f"{seq[pos : pos + line_width]}\n")
            else:
                writer.write(f"{seq}\n")

    if isinstance(output, str | Path):
        with open(output, "w", encoding="utf-8") as f:
            _write(f)
    else:
        _write(output)


def emit_excluded_fasta(
    graph: Any,
    resolutions: Sequence[Any],
    output: Path | TextIO | str,
    line_width: int | None = None,
    placed_unitigs: Iterable[str] | None = None,
) -> int:
    """Emit all unitigs excluded from scaffolds via excluded_unitig_ids to FASTA format.

    Headers include provenance tracking:
    >{unitig_id} excluded_by={resolution_feature_ids} length={length}

    Args:
        graph: Populated assembly graph containing segment sequences.
        resolutions: Sequence of ResolutionRecord objects.
        output: Destination file path or open text stream.
        line_width: Optional line wrap width (None for unwrapped).
        placed_unitigs: Optional iterable of unitig IDs that were placed in scaffolds.
            Any unitig placed in scaffolds is strictly excluded from excluded.fasta.

    Returns:
        Number of distinct excluded unitigs written.
    """
    placed_set = set(placed_unitigs) if placed_unitigs is not None else set()
    excluded_by: dict[str, list[str]] = defaultdict(list)
    for rec in resolutions:
        rec_fid = getattr(rec, "feature_id", str(rec))
        for uid in getattr(rec, "excluded_unitig_ids", ()):
            if uid not in placed_set:
                excluded_by[uid].append(rec_fid)

    sorted_uids = sorted(excluded_by.keys())

    def _write(writer: TextIO) -> None:
        for uid in sorted_uids:
            seq = ""
            if hasattr(graph, "get_sequence"):
                seq = graph.get_sequence(uid) or ""
            if not seq and hasattr(graph, "segments") and uid in graph.segments:
                seq = graph.segments[uid].sequence or ""
            seq_len = (
                len(seq)
                if seq
                else (graph.segments[uid].length if hasattr(graph, "segments") and uid in graph.segments else 0)
            )
            excl_str = ",".join(dict.fromkeys(excluded_by[uid]))
            header = f">{uid} excluded_by={excl_str} length={seq_len}"
            writer.write(f"{header}\n")
            if line_width is not None and line_width > 0:
                for pos in range(0, len(seq), line_width):
                    writer.write(f"{seq[pos : pos + line_width]}\n")
            else:
                writer.write(f"{seq}\n")

    if isinstance(output, str | Path):
        with open(output, "w", encoding="utf-8") as f:
            _write(f)
    else:
        _write(output)

    return len(sorted_uids)


def emit_agp(
    records: Sequence[AgpRecord],
    output: Path | TextIO | str,
) -> None:
    """Emit AGP v2.1 records with empirical provenance metadata tags.

    Args:
        records: Sequence of AgpRecord objects.
        output: Destination file path or open text stream.
    """

    def _write(writer: TextIO) -> None:
        writer.write("##agp-version 2.1\n")
        for rec in records:
            writer.write(rec.to_line() + "\n")

    if isinstance(output, str | Path):
        with open(output, "w", encoding="utf-8") as f:
            _write(f)
    else:
        _write(output)


def parse_agp(file: Path | TextIO | str) -> list[AgpRecord]:
    """Parse an AGP v2.1 file or stream containing provenance metadata tags.

    Args:
        file: Destination file path, open text stream, or AGP string content.

    Returns:
        List of deserialized AgpRecord objects.
    """

    def _parse(stream: Iterable[str]) -> list[AgpRecord]:
        records: list[AgpRecord] = []
        for raw_line in stream:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue

            comment_str: str | None = None
            prov: SequenceProvenance | None = None
            unpolished = False

            if "#" in line:
                data_part, comment_part = line.split("#", 1)
                data_tokens = [t.strip() for t in data_part.split("\t") if t.strip()]
                c_text = comment_part.strip()
                for chunk in c_text.split(";"):
                    chunk = chunk.strip()
                    if chunk.startswith("provenance="):
                        prov_val = chunk.split("=", 1)[1].strip()
                        prov = SequenceProvenance(prov_val)
                    elif chunk.startswith("unpolished="):
                        unpol_val = chunk.split("=", 1)[1].strip().lower()
                        unpolished = unpol_val in ("true", "1")
                    elif chunk:
                        comment_str = chunk if comment_str is None else f"{comment_str}; {chunk}"
            else:
                data_tokens = [t.strip() for t in line.split("\t") if t.strip()]

            (
                obj_id,
                obj_beg_s,
                obj_end_s,
                part_no_s,
                comp_type,
                comp_id,
                comp_beg_s,
                comp_end_s,
                orientation,
                *_,
            ) = data_tokens

            comp_beg: int | str
            comp_end: int | str
            try:
                comp_beg = int(comp_beg_s)
            except ValueError:
                comp_beg = comp_beg_s
            try:
                comp_end = int(comp_end_s)
            except ValueError:
                comp_end = comp_end_s

            records.append(
                AgpRecord(
                    object_id=obj_id,
                    object_beg=int(obj_beg_s),
                    object_end=int(obj_end_s),
                    part_number=int(part_no_s),
                    component_type=comp_type,
                    component_id=comp_id,
                    component_beg=comp_beg,
                    component_end=comp_end,
                    orientation=orientation,
                    sequence_provenance=prov,
                    unpolished=unpolished,
                    comment=comment_str,
                )
            )
        return records

    if isinstance(file, Path):
        with open(file, encoding="utf-8") as f:
            return _parse(f)
    if isinstance(file, str):
        if "\n" not in file and "\t" not in file:
            p = Path(file)
            if p.is_file():
                with open(p, encoding="utf-8") as f:
                    return _parse(f)
        return _parse(file.splitlines())
    return _parse(file)


DEFAULT_MIN_SUPER_LENGTH = int("100000")
PERCENT_MULTIPLIER = float("100.0")
TABLE_SEPARATOR_WIDTH = int("105")


def emit_grs_bed(
    agp_records: Sequence[AgpRecord | AgpGapRow],
    out_path: Path | str,
) -> Path:
    """Emit 0-based half-open BED records for eliminated germline-retained sequences (GRS).

    Strictly excludes variable repeats and satellite gaps (e.g. gap components 'N', 'U',
    gap_type == 'repeat', ONT-derived sequences, or SAT_GAP_VARIABLE).

    Args:
        agp_records: Assembled scaffold AGP v2.1 records.
        out_path: Destination BED file path.

    Returns:
        Path to the written BED file.
    """
    target = Path(out_path)
    target.parent.mkdir(parents=True, exist_ok=True)

    grs_provs = {
        SequenceProvenance.HIFI_GERMLINE,
        SequenceProvenance.hifi_germline,
        "HiFi-germline",
    }
    excluded_provs = {
        SequenceProvenance.ONT_DERIVED,
        SequenceProvenance.ont_derived,
        "ONT-derived",
    }

    bed_rows: list[tuple[str, int, int, str, int, str]] = []

    for rec in agp_records:
        comp_type = getattr(rec, "component_type", "")
        # Strictly exclude variable repeat gaps and non-sequence gap rows
        if comp_type in ("N", "U") or getattr(rec, "is_gap", False):
            continue
        if getattr(rec, "gap_type", None) == "repeat":
            continue

        prov = getattr(rec, "sequence_provenance", None)
        if prov in excluded_provs:
            continue
        if prov not in grs_provs:
            continue

        comment = getattr(rec, "comment", "") or ""
        if "SAT_GAP_VARIABLE" in comment or "repeat" in comment:
            continue

        chrom = rec.object_id
        chrom_start = rec.object_beg - 1
        chrom_end = rec.object_end
        name = str(rec.component_id)
        score = 0
        strand = rec.orientation if rec.orientation in ("+", "-") else "+"

        bed_rows.append((chrom, chrom_start, chrom_end, name, score, strand))

    bed_rows.sort(key=lambda r: (r[0], r[1], r[2]))

    with open(target, "w", encoding="utf-8") as f:
        for chrom, chrom_start, chrom_end, name, score, strand in bed_rows:
            f.write(f"{chrom}\t{chrom_start}\t{chrom_end}\t{name}\t{score}\t{strand}\n")

    return target


def emit_rdna_bed(
    rdna_regions: Sequence[Any],
    out_path: Path | str,
) -> Path:
    """Emit 0-based half-open BED records for defined 45S rRNA regions.

    Args:
        rdna_regions: Sequence of RdnaRegion instances or objects with
            sequence_id, start, end, hit_count, strand, and features/to_bed_line.
        out_path: Destination BED file path.

    Returns:
        Path to the written BED file.
    """
    target = Path(out_path)
    target.parent.mkdir(parents=True, exist_ok=True)

    sorted_regions = sorted(
        rdna_regions,
        key=lambda r: (getattr(r, "sequence_id", ""), getattr(r, "start", 0), getattr(r, "end", 0)),
    )
    with open(target, "w", encoding="utf-8") as f:
        for r in sorted_regions:
            if hasattr(r, "to_bed_line"):
                f.write(r.to_bed_line())
            else:
                seq_id = getattr(r, "sequence_id", "")
                s = getattr(r, "start", 0)
                e = getattr(r, "end", 0)
                st = getattr(r, "strand", Strand.FORWARD)
                st_sym = st.to_symbol() if hasattr(st, "to_symbol") else str(st)
                cnt = getattr(r, "hit_count", 0)
                f.write(f"{seq_id}\t{s}\t{e}\t45S_rRNA_region\t{cnt}\t{st_sym}\n")

    return target


def emit_telomere_bed(
    scaffolds: Mapping[str, str],
    out_path: Path | str,
    policy: Policy = DEFAULT_POLICY,
) -> Path:
    """Emit authentic germline terminal telomere cap coordinates for assembled scaffolds in BED format.

    Uses scan_terminal_tract to inspect the 5' prefix (testing motif_reverse, fallback motif_forward)
    and 3' suffix (testing motif_forward, fallback motif_reverse) while tolerating biological noise.

    Args:
        scaffolds: Mapping from scaffold/chromosome identifier to assembled nucleotide sequence.
        out_path: Destination BED file path.
        policy: Finish policy supplying canonical motifs and tolerance thresholds.

    Returns:
        Path to the emitted BED file.
    """
    target = Path(out_path)
    target.parent.mkdir(parents=True, exist_ok=True)

    tel_policy = policy.telomere
    motif_fwd = tel_policy.motif_forward.value
    motif_rev = tel_policy.motif_reverse.value
    terminal_scan_bp = tel_policy.terminal_scan_bp.value
    max_offset = tel_policy.max_terminal_offset_bp.value
    max_internal_gap = tel_policy.max_internal_gap_bp.value
    min_repeats = tel_policy.min_terminal_repeat_count.value

    bed_rows: list[tuple[str, int, int, str, int, str]] = []

    for scaf_id, seq in scaffolds.items():
        if not seq:
            continue
        seq_len = len(seq)

        # 5' prefix: scan motif_reverse first, fallback to motif_forward
        call_5p = scan_terminal_tract(
            seq,
            motif=motif_rev,
            terminal_scan_bp=terminal_scan_bp,
            max_terminal_offset_bp=max_offset,
            max_internal_gap_bp=max_internal_gap,
            min_repeat_count=min_repeats,
            is_suffix=False,
        )
        matched_5p_motif = motif_rev
        if not call_5p.is_terminal_cap:
            call_5p_fwd = scan_terminal_tract(
                seq,
                motif=motif_fwd,
                terminal_scan_bp=terminal_scan_bp,
                max_terminal_offset_bp=max_offset,
                max_internal_gap_bp=max_internal_gap,
                min_repeat_count=min_repeats,
                is_suffix=False,
            )
            if call_5p_fwd.is_terminal_cap:
                call_5p = call_5p_fwd
                matched_5p_motif = motif_fwd

        if call_5p.is_terminal_cap:
            chrom_start = call_5p.offset_bp
            chrom_end = call_5p.offset_bp + call_5p.tract_length_bp
            strand = "-" if matched_5p_motif == motif_rev else "+"
            bed_rows.append((scaf_id, chrom_start, chrom_end, "telomere_5p", call_5p.repeat_count, strand))

        # 3' suffix: scan motif_forward first, fallback to motif_reverse
        call_3p = scan_terminal_tract(
            seq,
            motif=motif_fwd,
            terminal_scan_bp=terminal_scan_bp,
            max_terminal_offset_bp=max_offset,
            max_internal_gap_bp=max_internal_gap,
            min_repeat_count=min_repeats,
            is_suffix=True,
        )
        matched_3p_motif = motif_fwd
        if not call_3p.is_terminal_cap:
            call_3p_rev = scan_terminal_tract(
                seq,
                motif=motif_rev,
                terminal_scan_bp=terminal_scan_bp,
                max_terminal_offset_bp=max_offset,
                max_internal_gap_bp=max_internal_gap,
                min_repeat_count=min_repeats,
                is_suffix=True,
            )
            if call_3p_rev.is_terminal_cap:
                call_3p = call_3p_rev
                matched_3p_motif = motif_rev

        if call_3p.is_terminal_cap:
            chrom_end = seq_len - call_3p.offset_bp
            chrom_start = chrom_end - call_3p.tract_length_bp
            strand = "+" if matched_3p_motif == motif_fwd else "-"
            bed_rows.append((scaf_id, chrom_start, chrom_end, "telomere_3p", call_3p.repeat_count, strand))

    bed_rows.sort(key=lambda r: (r[0], r[1], r[2]))

    with open(target, "w", encoding="utf-8") as f:
        for chrom, chrom_start, chrom_end, name, score, strand in bed_rows:
            f.write(f"{chrom}\t{chrom_start}\t{chrom_end}\t{name}\t{score}\t{strand}\n")

    return target


def find_break_site_scaffolds(
    scaffolds: Mapping[str, str],
    agp_records: Sequence[AgpRecord | AgpGapRow],
    junctions: Sequence[Any],
) -> dict[str, list[Any]]:
    """Map each scaffold ID to its matching break site junction records.

    Matches either j.stem_port, j.germline_continuation_port, or j.chrom
    against the scaffold's constituent segment IDs.
    """
    scaf_segments: dict[str, set[str]] = defaultdict(set)
    for r in agp_records:
        if getattr(r, "component_type", "") == "W":
            scaf_segments[r.object_id].add(str(r.component_id))

    matches_by_scaf: dict[str, list[Any]] = defaultdict(list)
    for scaf_id in scaffolds:
        segs = scaf_segments.get(scaf_id, set()) | {scaf_id}
        for j in junctions:
            j_stem = getattr(j, "stem_port", "") or ""
            j_stem_clean = j_stem.rstrip("+-")
            j_germ = getattr(j, "germline_continuation_port", "") or ""
            j_germ_clean = j_germ.rstrip("+-") if j_germ else ""
            j_chrom = getattr(j, "chrom", "") or ""

            if (
                j_chrom in segs
                or j_stem_clean in segs
                or j_stem in segs
                or (j_germ_clean and j_germ_clean in segs)
                or (j_germ and j_germ in segs)
            ):
                matches_by_scaf[scaf_id].append(j)

    return matches_by_scaf


def nominate_super_scaffolds(
    scaffolds: Mapping[str, str],
    agp_records: Sequence[AgpRecord | AgpGapRow],
    junctions: Sequence[Any],
    min_super_length: int = DEFAULT_MIN_SUPER_LENGTH,
    resolutions: list[Any] | None = None,
    end_states: list[Any] | None = None,
) -> tuple[dict[str, str], list[AgpRecord | AgpGapRow], dict[str, str]]:
    """Nominate large scaffolds carrying break site evidence as candidate chromosomes (SUPER_N).

    Candidates are ranked in descending order of sequence length: SUPER_1 for largest,
    SUPER_2 for second largest, etc. Non-candidate sequences retain their original IDs.

    Args:
        scaffolds: Mapping of scaffold IDs to nucleotide sequences.
        agp_records: Assembled AGP records.
        junctions: Sequence of JunctionRecord objects representing break sites.
        min_super_length: Minimum sequence length in bp to qualify as candidate chromosome.
        resolutions: Optional list of ResolutionRecord objects to update in-place.
        end_states: Optional list of EndStateRecord objects to update in-place.

    Returns:
        Tuple of (renamed_scaffolds, updated_agp_records, rename_map).
    """
    matches = find_break_site_scaffolds(scaffolds, agp_records, junctions)
    qualifying = [
        sid for sid, seq in scaffolds.items() if len(matches.get(sid, [])) >= 1 and len(seq) >= min_super_length
    ]
    sorted_candidates = sorted(qualifying, key=lambda sid: (-len(scaffolds[sid]), sid))

    rename_map: dict[str, str] = {}
    for idx, sid in enumerate(sorted_candidates, start=1):
        rename_map[sid] = f"SUPER_{idx}"

    new_scaffolds: dict[str, str] = {}
    # Candidate chromosomes first, sorted
    for sid in sorted_candidates:
        new_scaffolds[rename_map[sid]] = scaffolds[sid]
    # Remaining non-candidate scaffolds
    for sid, seq in scaffolds.items():
        if sid not in rename_map:
            new_scaffolds[sid] = seq

    new_agp: list[AgpRecord | AgpGapRow] = []
    for rec in agp_records:
        if rec.object_id in rename_map:
            rec.object_id = rename_map[rec.object_id]
        new_agp.append(rec)

    if resolutions is not None:
        for rec in resolutions:
            if getattr(rec, "chrom", None) in rename_map:
                rec.chrom = rename_map[rec.chrom]

    if end_states is not None:
        for rec in end_states:
            if getattr(rec, "chrom", None) in rename_map:
                rec.chrom = rename_map[rec.chrom]

    return new_scaffolds, new_agp, rename_map


@dataclass(slots=True)
class ResolutionSummaryMetrics:
    """Tally of resolved complexities / quirks across the assembly."""

    somatic_allelic_bubbles: int = 0
    variable_satellite_gaps: int = 0
    fixed_satellite_gaps: int = 0
    intercalary_grs_bridges: int = 0
    pde_terminal_caps: int = 0
    subtelomeric_satellite_caps: int = 0
    rdna_cassettes_traversed: int = 0
    heuristic_caps_closed: int = 0


@dataclass(slots=True)
class SequenceSummaryMetrics:
    """Stoichiometric and sequence metrics for a break-site qualified chromosome/scaffold."""

    sequence_id: str
    core_length_bp: int
    core_depth: float
    grs_length_bp: int
    grs_depth: float
    stoichiometry_ratio: float
    grs_percentage: float
    break_site_count: int
    somatic_allelic_bubbles: int = 0
    variable_satellite_gaps: int = 0
    fixed_satellite_gaps: int = 0
    intercalary_grs_bridges: int = 0
    pde_terminal_caps: int = 0
    subtelomeric_satellite_caps: int = 0
    rdna_cassettes_traversed: int = 0
    heuristic_caps_closed: int = 0


@dataclass(slots=True)
class BreakSiteGenomeSummary:
    """Summary metrics across all break-site verified candidate chromosomes."""

    qualifying_sequences: list[SequenceSummaryMetrics]
    total_core_length_bp: int
    mean_core_depth: float
    total_grs_length_bp: int
    mean_grs_depth: float
    overall_stoichiometry_ratio: float
    overall_grs_percentage: float
    total_break_sites: int
    resolution_summary: ResolutionSummaryMetrics = field(default_factory=ResolutionSummaryMetrics)


def compute_break_site_summary(
    graph: Any,
    scaffolds: Mapping[str, str],
    agp_records: Sequence[AgpRecord | AgpGapRow],
    junctions: Sequence[Any],
    resolutions: Sequence[Any] = (),
) -> BreakSiteGenomeSummary:
    """Compute core genome vs GRS coverage strictly across sequences with break site evidence.

    Args:
        graph: Bidirected assembly graph containing segment depths.
        scaffolds: Mapping of scaffold ID to nucleotide sequence.
        agp_records: Assembled AGP records.
        junctions: Junction records representing break sites.
        resolutions: Optional sequence of ResolutionRecord objects to tally resolved quirks.

    Returns:
        BreakSiteGenomeSummary with per-sequence breakdown and genome-wide totals.
    """
    matches = find_break_site_scaffolds(scaffolds, agp_records, junctions)
    qualifying_ids = [sid for sid in scaffolds if len(matches.get(sid, [])) >= 1]
    sorted_qualifying = sorted(qualifying_ids, key=lambda sid: (-len(scaffolds[sid]), sid))

    agp_by_scaf: dict[str, list[AgpRecord | AgpGapRow]] = defaultdict(list)
    scaf_segs: dict[str, set[str]] = defaultdict(set)
    for rec in agp_records:
        agp_by_scaf[rec.object_id].append(rec)
        if getattr(rec, "component_type", "") == "W":
            scaf_segs[rec.object_id].add(str(rec.component_id))

    core_provs = {
        SequenceProvenance.HIFI_CORE,
        SequenceProvenance.hifi_core,
        "HiFi-core",
    }
    grs_provs = {
        SequenceProvenance.HIFI_GERMLINE,
        SequenceProvenance.hifi_germline,
        "HiFi-germline",
    }

    # Map resolutions to scaffolds
    res_by_scaf: dict[str, list[Any]] = defaultdict(list)
    for r in resolutions:
        r_chrom = getattr(r, "chrom", None) or ""
        matched = False
        for sid in scaffolds:
            if r_chrom == sid or r_chrom in scaf_segs.get(sid, set()):
                res_by_scaf[sid].append(r)
                matched = True
                break
        if not matched and r_chrom:
            res_by_scaf[r_chrom].append(r)

    def _tally_resolutions(res_list: Sequence[Any]) -> ResolutionSummaryMetrics:
        tally = ResolutionSummaryMetrics()
        for rec in res_list:
            rtype = str(getattr(rec, "resolution_type", "") or "")
            if rtype in ("somatic_allelic_bubble", "ResolutionType.SOMATIC_ALLELIC_BUBBLE"):
                tally.somatic_allelic_bubbles += 1
            elif rtype in ("sat_gap_variable", "ResolutionType.SAT_GAP_VARIABLE"):
                tally.variable_satellite_gaps += 1
            elif rtype in ("sat_gap_fixed", "ResolutionType.SAT_GAP_FIXED"):
                tally.fixed_satellite_gaps += 1
            elif rtype in ("intercalary_grs_bridge", "ResolutionType.INTERCALARY_GRS_BRIDGE"):
                tally.intercalary_grs_bridges += 1
            elif rtype in ("pde_terminal_cap", "ResolutionType.PDE_TERMINAL_CAP"):
                tally.pde_terminal_caps += 1
            elif rtype in ("subtelomeric_satellite_cap", "ResolutionType.SUBTELOMERIC_SATELLITE_CAP"):
                tally.subtelomeric_satellite_caps += 1
            elif rtype in ("rdna_cassette", "ResolutionType.RDNA_CASSETTE"):
                tally.rdna_cassettes_traversed += 1

            if getattr(rec, "capped_via_heuristic", False) or "heuristic" in str(getattr(rec, "notes", "")).lower():
                tally.heuristic_caps_closed += 1
        return tally

    qualifying_metrics: list[SequenceSummaryMetrics] = []
    total_core_len_all = 0
    total_core_weighted_depth_all = 0.0
    total_grs_len_all = 0
    total_grs_weighted_depth_all = 0.0
    total_breaks_all = 0

    for sid in sorted_qualifying:
        scaf_recs = agp_by_scaf.get(sid, [])
        core_len = 0
        core_weighted_depth = 0.0
        grs_len = 0
        grs_weighted_depth = 0.0

        for r in scaf_recs:
            if getattr(r, "component_type", "") != "W":
                continue
            comp_len = r.object_end - r.object_beg + 1
            seg_id = str(r.component_id)
            seg_depth = 0.0
            if hasattr(graph, "segments") and seg_id in graph.segments:
                d = graph.segments[seg_id].depth
                seg_depth = float(d) if d is not None else 0.0

            prov = getattr(r, "sequence_provenance", None)
            if prov in core_provs:
                core_len += comp_len
                core_weighted_depth += comp_len * seg_depth
            elif prov in grs_provs:
                grs_len += comp_len
                grs_weighted_depth += comp_len * seg_depth

        core_depth = (core_weighted_depth / core_len) if core_len > 0 else 0.0
        grs_depth = (grs_weighted_depth / grs_len) if grs_len > 0 else 0.0
        stoich_ratio = (grs_depth / core_depth) if core_depth > 0.0 else 0.0
        denom = core_len + grs_len
        grs_pct = (grs_len / denom * PERCENT_MULTIPLIER) if denom > 0 else 0.0
        n_breaks = len(matches.get(sid, []))
        scaf_tally = _tally_resolutions(res_by_scaf.get(sid, []))

        qualifying_metrics.append(
            SequenceSummaryMetrics(
                sequence_id=sid,
                core_length_bp=core_len,
                core_depth=core_depth,
                grs_length_bp=grs_len,
                grs_depth=grs_depth,
                stoichiometry_ratio=stoich_ratio,
                grs_percentage=grs_pct,
                break_site_count=n_breaks,
                somatic_allelic_bubbles=scaf_tally.somatic_allelic_bubbles,
                variable_satellite_gaps=scaf_tally.variable_satellite_gaps,
                fixed_satellite_gaps=scaf_tally.fixed_satellite_gaps,
                intercalary_grs_bridges=scaf_tally.intercalary_grs_bridges,
                pde_terminal_caps=scaf_tally.pde_terminal_caps,
                subtelomeric_satellite_caps=scaf_tally.subtelomeric_satellite_caps,
                rdna_cassettes_traversed=scaf_tally.rdna_cassettes_traversed,
                heuristic_caps_closed=scaf_tally.heuristic_caps_closed,
            )
        )

        total_core_len_all += core_len
        total_core_weighted_depth_all += core_weighted_depth
        total_grs_len_all += grs_len
        total_grs_weighted_depth_all += grs_weighted_depth
        total_breaks_all += n_breaks

    overall_core_depth = (total_core_weighted_depth_all / total_core_len_all) if total_core_len_all > 0 else 0.0
    overall_grs_depth = (total_grs_weighted_depth_all / total_grs_len_all) if total_grs_len_all > 0 else 0.0
    overall_stoich_ratio = (overall_grs_depth / overall_core_depth) if overall_core_depth > 0.0 else 0.0
    overall_denom = total_core_len_all + total_grs_len_all
    overall_grs_pct = (total_grs_len_all / overall_denom * PERCENT_MULTIPLIER) if overall_denom > 0 else 0.0
    overall_tally = _tally_resolutions(resolutions)

    return BreakSiteGenomeSummary(
        qualifying_sequences=qualifying_metrics,
        total_core_length_bp=total_core_len_all,
        mean_core_depth=overall_core_depth,
        total_grs_length_bp=total_grs_len_all,
        mean_grs_depth=overall_grs_depth,
        overall_stoichiometry_ratio=overall_stoich_ratio,
        overall_grs_percentage=overall_grs_pct,
        total_break_sites=total_breaks_all,
        resolution_summary=overall_tally,
    )


def emit_summary_tsv(
    summary: BreakSiteGenomeSummary,
    out_path: Path | str,
) -> Path:
    """Emit core vs GRS summary metrics table to TSV format with resolution event tallies.

    Args:
        summary: Computed BreakSiteGenomeSummary.
        out_path: Destination TSV file path.

    Returns:
        Path to the emitted TSV file.
    """
    target = Path(out_path)
    target.parent.mkdir(parents=True, exist_ok=True)

    header = [
        "sequence_id",
        "core_length_bp",
        "core_depth",
        "grs_length_bp",
        "grs_depth",
        "stoichiometry_ratio",
        "grs_percentage",
        "break_site_count",
        "somatic_allelic_bubbles",
        "variable_satellite_gaps",
        "fixed_satellite_gaps",
        "intercalary_grs_bridges",
        "pde_terminal_caps",
        "subtelomeric_satellite_caps",
        "rdna_cassettes_traversed",
        "heuristic_caps_closed",
    ]

    lines = ["\t".join(header)]
    for seq in summary.qualifying_sequences:
        lines.append(
            f"{seq.sequence_id}\t{seq.core_length_bp}\t{seq.core_depth:.2f}\t"
            f"{seq.grs_length_bp}\t{seq.grs_depth:.2f}\t{seq.stoichiometry_ratio:.4f}\t"
            f"{seq.grs_percentage:.2f}\t{seq.break_site_count}\t"
            f"{seq.somatic_allelic_bubbles}\t{seq.variable_satellite_gaps}\t"
            f"{seq.fixed_satellite_gaps}\t{seq.intercalary_grs_bridges}\t"
            f"{seq.pde_terminal_caps}\t{seq.subtelomeric_satellite_caps}\t"
            f"{seq.rdna_cassettes_traversed}\t{seq.heuristic_caps_closed}"
        )

    # Summary TOTAL row
    res = summary.resolution_summary
    lines.append(
        f"TOTAL\t{summary.total_core_length_bp}\t{summary.mean_core_depth:.2f}\t"
        f"{summary.total_grs_length_bp}\t{summary.mean_grs_depth:.2f}\t"
        f"{summary.overall_stoichiometry_ratio:.4f}\t{summary.overall_grs_percentage:.2f}\t"
        f"{summary.total_break_sites}\t"
        f"{res.somatic_allelic_bubbles}\t{res.variable_satellite_gaps}\t"
        f"{res.fixed_satellite_gaps}\t{res.intercalary_grs_bridges}\t"
        f"{res.pde_terminal_caps}\t{res.subtelomeric_satellite_caps}\t"
        f"{res.rdna_cassettes_traversed}\t{res.heuristic_caps_closed}"
    )

    with open(target, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")

    return target


def format_summary_table(summary: BreakSiteGenomeSummary) -> str:
    """Format human-readable summary table for console output."""
    sep = "=" * TABLE_SEPARATOR_WIDTH
    sub_sep = "-" * TABLE_SEPARATOR_WIDTH
    header_title = "gfapdefinisher Finish Summary (Break-Site Verified Sequences)"
    table_lines = [
        sep,
        header_title,
        sep,
        f"{'Sequence ID':<18}{'Core (bp)':>14}{'Core Depth':>14}{'GRS (bp)':>14}{'GRS Depth':>14}{'Stoich Ratio':>16}{'GRS %':>10}{'Breaks':>8}",
        sub_sep,
    ]

    for s in summary.qualifying_sequences:
        table_lines.append(
            f"{s.sequence_id:<18}{s.core_length_bp:>14,}{s.core_depth:>13.2f}x{s.grs_length_bp:>14,}{s.grs_depth:>13.2f}x{s.stoichiometry_ratio:>16.4f}{s.grs_percentage:>9.2f}%{s.break_site_count:>8}"
        )

    table_lines.append(sub_sep)
    table_lines.append(
        f"{'TOTAL':<18}{summary.total_core_length_bp:>14,}{summary.mean_core_depth:>13.2f}x{summary.total_grs_length_bp:>14,}{summary.mean_grs_depth:>13.2f}x{summary.overall_stoichiometry_ratio:>16.4f}{summary.overall_grs_percentage:>9.2f}%{summary.total_break_sites:>8}"
    )
    table_lines.append(sep)

    res = summary.resolution_summary
    table_lines.append("Resolved Complexities (Quirks):")
    table_lines.append(f"  - Somatic allelic bubbles popped:        {res.somatic_allelic_bubbles:>5}")
    table_lines.append(f"  - Variable-size satellite gaps:          {res.variable_satellite_gaps:>5}")
    table_lines.append(f"  - Fixed-size satellite gaps:             {res.fixed_satellite_gaps:>5}")
    table_lines.append(f"  - Intercalary GRS bridges:               {res.intercalary_grs_bridges:>5}")
    table_lines.append(f"  - PDE terminal caps stitched:            {res.pde_terminal_caps:>5}")
    table_lines.append(f"  - Subtelomeric satellite caps:           {res.subtelomeric_satellite_caps:>5}")
    table_lines.append(f"  - rDNA cassettes traversed:              {res.rdna_cassettes_traversed:>5}")
    table_lines.append(f"  - Heuristic repeat port closures:        {res.heuristic_caps_closed:>5}")
    table_lines.append(sep)

    return "\n".join(table_lines)


def parse_summary_tsv(file: Path | str) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Parse a summary.tsv file into per-sequence records and the overall total record."""
    p = Path(file)
    content = p.read_text(encoding="utf-8")
    lines = [line.strip() for line in content.splitlines() if line.strip()]
    if not lines:
        return [], {}

    header = lines[0].split("\t")
    seq_rows: list[dict[str, Any]] = []
    total_row: dict[str, Any] = {}

    for line in lines[1:]:
        parts = line.split("\t")
        row = dict(zip(header, parts, strict=False))
        typed_row: dict[str, Any] = {
            "sequence_id": row["sequence_id"],
            "core_length_bp": int(row["core_length_bp"]),
            "core_depth": float(row["core_depth"]),
            "grs_length_bp": int(row["grs_length_bp"]),
            "grs_depth": float(row["grs_depth"]),
            "stoichiometry_ratio": float(row["stoichiometry_ratio"]),
            "grs_percentage": float(row["grs_percentage"]),
            "break_site_count": int(row["break_site_count"]),
            "somatic_allelic_bubbles": int(row.get("somatic_allelic_bubbles", 0)),
            "variable_satellite_gaps": int(row.get("variable_satellite_gaps", 0)),
            "fixed_satellite_gaps": int(row.get("fixed_satellite_gaps", 0)),
            "intercalary_grs_bridges": int(row.get("intercalary_grs_bridges", 0)),
            "pde_terminal_caps": int(row.get("pde_terminal_caps", 0)),
            "subtelomeric_satellite_caps": int(row.get("subtelomeric_satellite_caps", 0)),
            "rdna_cassettes_traversed": int(row.get("rdna_cassettes_traversed", 0)),
            "heuristic_caps_closed": int(row.get("heuristic_caps_closed", 0)),
        }
        if typed_row["sequence_id"] == "TOTAL":
            total_row = typed_row
        else:
            seq_rows.append(typed_row)

    return seq_rows, total_row
