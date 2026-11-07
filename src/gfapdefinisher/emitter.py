"""Sequence and AGP assembly emitter with empirical provenance tracking.

Stitches unitig paths, rDNA cassettes, and bridged gaps into finished FASTA
and AGP v2.1 outputs with empirical provenance metadata and unpolished sequence flags.
Implements BIO-004 (repeat preservation) and BIO-005 (lossless overlap trimming).
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any, TextIO

from gfapdefinisher.models import Strand
from gfapdefinisher.policy import DEFAULT_POLICY, Policy
from gfapdefinisher.policy.models import PolicyValue

__all__ = [
    "AgpGapRow",
    "AgpPolicy",
    "AgpRecord",
    "DEFAULT_AGP_POLICY",
    "FastaAgpEmitter",
    "PathComponent",
    "ScaffoldPlan",
    "ScaffoldStep",
    "SequenceProvenance",
    "assemble_path_sequence",
    "assemble_scaffold_sequence",
    "emit_agp",
    "emit_fasta",
    "emit_scaffold",
    "parse_agp",
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
        linkage_evidence: AGP linkage evidence (e.g. 'rna_tangle', 'paired-ends').
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
    linkage_evidence: str = "paired-ends"
    gap_component_type: str = "N"

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
        linkage_evidence: str = "paired-ends",
        gap_component_type: str = "N",
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
                and comp.linkage_evidence in ("rna_tangle", "paired-ends")
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
