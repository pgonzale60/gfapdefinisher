"""Biological quirk, junction, and chromosome end-state taxonomy.

Defines the ratified 3-layer biological feature schema for assembly graph
finishing, capturing chromosome end completeness, somatic diminution (PDE)
junctions, and interval-level structural quirks.

Ratified in Issue #26.
"""

from __future__ import annotations

import csv
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import asdict, dataclass, field
from enum import Enum
from io import StringIO
from pathlib import Path
from typing import Any, TextIO

from gfapdefinisher.classifier import (
    PdeJunctionCall,
    _lookup_evidence,
    classify_pde_junction,
)
from gfapdefinisher.graph import BidirectedGraph
from gfapdefinisher.ledger import EndEvidence
from gfapdefinisher.models import SegmentPort, Strand
from gfapdefinisher.policy import DEFAULT_POLICY, Policy

__all__ = [
    "END_STATE_TSV_COLUMNS",
    "JUNCTION_TSV_COLUMNS",
    "QUIRK_COLUMNS",
    "QUIRK_TSV_COLUMNS",
    "RESOLUTION_COLUMNS",
    "RESOLUTION_TSV_COLUMNS",
    "EndState",
    "EndStateRecord",
    "IntervalQuirk",
    "JunctionClass",
    "JunctionRecord",
    "QuirkRecord",
    "ResolutionRecord",
    "ResolutionType",
    "build_feature_id",
    "collect_features",
    "create_junction_record",
    "emit_end_states_tsv",
    "emit_junctions_tsv",
    "emit_quirks_tsv",
    "emit_resolutions_tsv",
    "parse_end_states_tsv",
    "parse_junctions_tsv",
    "parse_quirks_tsv",
    "parse_resolutions_tsv",
]


class EndState(str, Enum):
    """Chromosome end-state classification across somatic and germline lineages."""

    END_COMPLETE = "end_complete"
    END_SOMATIC_MISSING = "end_somatic_missing"
    END_GERMLINE_ORPHAN = "end_germline_orphan"
    END_SHARED_EXEMPT = "end_shared_exempt"

    # Lowercase aliases for direct attribute access
    end_complete = "end_complete"
    end_somatic_missing = "end_somatic_missing"
    end_germline_orphan = "end_germline_orphan"
    end_shared_exempt = "end_shared_exempt"


class JunctionClass(str, Enum):
    """Classification of Programmed DNA Elimination (PDE) diminution junctions."""

    PDE_TERMINAL = "pde_terminal"
    PDE_INTERCALARY = "pde_intercalary"

    # Lowercase aliases
    pde_terminal = "pde_terminal"
    pde_intercalary = "pde_intercalary"


class IntervalQuirk(str, Enum):
    """Interval-level biological quirk or assembly graph idiosyncrasy."""

    SAT_GAP_VARIABLE = "sat_gap_variable"
    SAT_GAP_FIXED = "sat_gap_fixed"
    RDNA_CASSETTE = "rdna_cassette"
    UNCOLLAPSED_COPY = "uncollapsed_copy"
    LOW_SUPPORT_STUB = "low_support_stub"
    BUBBLE_REJOINING = "bubble_rejoining"
    ITR = "itr"
    SUBTELOMERIC_DUP = "subtelomeric_dup"
    INTERCALARY_ELIMINATION_CANDIDATE = "intercalary_elimination_candidate"
    SOMATIC_ALLELIC_BUBBLE = "somatic_allelic_bubble"
    INTERCALARY_GRS_BRIDGE = "intercalary_grs_bridge"
    PDE_TERMINAL_CAP = "pde_terminal_cap"
    SUBTELOMERIC_SATELLITE_CAP = "subtelomeric_satellite_cap"

    # Lowercase aliases
    sat_gap_variable = "sat_gap_variable"
    sat_gap_fixed = "sat_gap_fixed"
    rdna_cassette = "rdna_cassette"
    uncollapsed_copy = "uncollapsed_copy"
    low_support_stub = "low_support_stub"
    bubble_rejoining = "bubble_rejoining"
    itr = "itr"
    subtelomeric_dup = "subtelomeric_dup"
    intercalary_elimination_candidate = "intercalary_elimination_candidate"
    somatic_allelic_bubble = "somatic_allelic_bubble"
    intercalary_grs_bridge = "intercalary_grs_bridge"
    pde_terminal_cap = "pde_terminal_cap"
    subtelomeric_satellite_cap = "subtelomeric_satellite_cap"


def build_feature_id(
    feature_class: str | Enum,
    chrom: str,
    side: str,
    index: int | None = None,
) -> str:
    """Format a canonical feature identifier: `<feature_class>:<chrom>:<side>[:<index>]`.

    Args:
        feature_class: Class enum or string identifier.
        chrom: Chromosome name (e.g. 'I', 'IV', 'X').
        side: Arm end or locus designation (e.g. '5p', '3p', 'mid').
        index: Optional 1-based occurrence index.

    Returns:
        Formatted canonical identifier string.
    """
    cls_name = feature_class.value if isinstance(feature_class, Enum) else str(feature_class)
    if index is None:
        return f"{cls_name}:{chrom}:{side}"
    return f"{cls_name}:{chrom}:{side}:{index}"


@dataclass(slots=True)
class EndStateRecord:
    """Layer 1 record evaluating somatic and germline completeness for one chromosome end.

    Attributes:
        feature_id: Canonical identifier (e.g. 'end_complete:I:5p').
        chrom: Chromosome name.
        side: Chromosome arm end ('5p' or '3p').
        state: EndState enum classification.
        somatic_port: Oriented unitig port capping somatic lineage, if present.
        germline_port: Oriented unitig port capping germline lineage, if present.
        shared_port: Oriented unitig port when end is exempt from diminution.
        core_coverage_bp_median: Median read depth of core somatic segments.
        somatic_depth_ratio: Normalized somatic cap depth ratio relative to core.
        germline_depth_ratio: Normalized germline cap depth ratio relative to core.
        resolved: Whether both lineages are satisfactorily terminated.
        notes: Contextual provenance notes.
    """

    feature_id: str
    chrom: str
    side: str
    state: EndState
    somatic_port: str | None = None
    germline_port: str | None = None
    shared_port: str | None = None
    core_coverage_bp_median: float | None = None
    somatic_depth_ratio: float | None = None
    germline_depth_ratio: float | None = None
    resolved: bool = True
    notes: str = ""

    def __post_init__(self) -> None:
        if not self.chrom:
            raise ValueError("Chromosome cannot be empty.")
        if not self.side:
            raise ValueError("Side cannot be empty.")
        if isinstance(self.state, str) and not isinstance(self.state, EndState):
            self.state = EndState(self.state)

    def to_dict(self) -> dict[str, Any]:
        """Serialize record to dictionary with string enum representations."""
        d = asdict(self)
        d["state"] = self.state.value
        return d

    def to_tsv_row(self) -> list[str]:
        """Serialize record to an ordered list of string values for TSV emission."""
        return [
            self.feature_id,
            self.chrom,
            self.side,
            self.state.value,
            self.somatic_port or "",
            self.germline_port or "",
            self.shared_port or "",
            str(self.core_coverage_bp_median) if self.core_coverage_bp_median is not None else "",
            str(self.somatic_depth_ratio) if self.somatic_depth_ratio is not None else "",
            str(self.germline_depth_ratio) if self.germline_depth_ratio is not None else "",
            str(self.resolved),
            self.notes,
        ]


@dataclass(slots=True)
class JunctionRecord:
    """Layer 2 record classifying a Programmed DNA Elimination (PDE) fork or break.

    Attributes:
        feature_id: Canonical identifier (e.g. 'pde_terminal:IV:3p', 'pde_intercalary:X:mid:1').
        chrom: Chromosome name.
        side: Arm end or locus ('5p', '3p', 'mid').
        junction_class: JunctionClass enum classification.
        stem_port: Oriented unitig port entering the diminution boundary.
        somatic_terminus_port: Somatic branch port terminating or stubbing out.
        germline_continuation_port: Germline branch port continuing chromosome arm.
        break_index: 1-based index for intercalary breaks.
        core_depth_ratio: Depth ratio of germline branch relative to core somatic coverage ($C_{\text{core}}$).
        read_count: Supporting read count across the junction boundary.
        confidence_tier: Confidence tier designation (e.g. 'CONFIDENT').
        notes: Contextual provenance notes.
    """

    feature_id: str
    chrom: str
    side: str
    junction_class: JunctionClass
    stem_port: str
    somatic_terminus_port: str | None = None
    germline_continuation_port: str | None = None
    break_index: int = 1
    core_depth_ratio: float | None = None
    read_count: int | None = None
    confidence_tier: str = "CONFIDENT"
    notes: str = ""

    def __post_init__(self) -> None:
        if not self.chrom:
            raise ValueError("Chromosome cannot be empty.")
        if not self.side:
            raise ValueError("Side cannot be empty.")
        if not self.stem_port:
            raise ValueError("Stem port cannot be empty.")
        if self.break_index < 1:
            raise ValueError("Break index must be at least 1.")
        if isinstance(self.junction_class, str) and not isinstance(self.junction_class, JunctionClass):
            self.junction_class = JunctionClass(self.junction_class)

    def to_dict(self) -> dict[str, Any]:
        """Serialize record to dictionary with string enum representations."""
        d = asdict(self)
        d["junction_class"] = self.junction_class.value
        return d

    def to_tsv_row(self) -> list[str]:
        """Serialize record to an ordered list of string values for TSV emission."""
        return [
            self.feature_id,
            self.chrom,
            self.side,
            self.junction_class.value,
            self.stem_port,
            self.somatic_terminus_port or "",
            self.germline_continuation_port or "",
            str(self.break_index),
            str(self.core_depth_ratio) if self.core_depth_ratio is not None else "",
            str(self.read_count) if self.read_count is not None else "",
            self.confidence_tier,
            self.notes,
        ]


@dataclass(slots=True)
class QuirkRecord:
    """Layer 3 record documenting an interval-level biological quirk or assembly idiosyncrasy.

    Attributes:
        feature_id: Canonical identifier (e.g. 'sat_gap_variable:IV:3p:1').
        chrom: Chromosome name.
        side: Arm end or locus.
        quirk_type: IntervalQuirk enum classification.
        spanning_read_count: Number of reads anchoring both flanking unitigs.
        length_min_bp: Minimum observed interval length across spanning reads.
        length_median_bp: Median observed interval length across spanning reads.
        length_max_bp: Maximum observed interval length across spanning reads.
        length_iqr_bp: Interquartile range of observed interval length.
        repeat_monomer_len_bp: Length in bp of the repeating unit motif.
        copy_number_range: Empirical copy-number range string (e.g. '12-45').
        left_anchor_port: Bounding unitig port on the proximal side.
        right_anchor_port: Bounding unitig port on the distal side.
        agp_component_type: AGP component designation (e.g. 'W' or 'U').
        agp_linkage_evidence: AGP linkage evidence designation (e.g. 'ont-spanning').
        unitig_ids: Tuple of involved assembly graph unitig identifiers.
        read_count: Total distinct reads supporting this segment/feature.
        core_coverage_bp_median: Somatic core reference depth ($C_{core}$).
        core_depth_ratio: Feature depth relative to somatic core coverage ($D / C_{core}$).
        representative_sequence_selected: Whether a real read sequence was chosen (BIO-004).
        notes: Contextual notes or provenance information.
    """

    feature_id: str
    chrom: str
    side: str
    quirk_type: IntervalQuirk
    spanning_read_count: int | None = None
    length_min_bp: int | None = None
    length_median_bp: int | None = None
    length_max_bp: int | None = None
    length_iqr_bp: int | None = None
    repeat_monomer_len_bp: int | None = None
    copy_number_range: str | None = None
    left_anchor_port: str | None = None
    right_anchor_port: str | None = None
    agp_component_type: str | None = None
    agp_linkage_evidence: str | None = None
    unitig_ids: tuple[str, ...] = field(default_factory=tuple)
    read_count: int | None = None
    core_coverage_bp_median: float | None = None
    core_depth_ratio: float | None = None
    representative_sequence_selected: bool = False
    sequence_provenance: str | None = None
    unpolished: bool = False
    notes: str = ""
    agp_start: int | None = None
    agp_end: int | None = None
    excluded_unitig_ids: tuple[str, ...] = field(default_factory=tuple)

    def __post_init__(self) -> None:
        if not self.chrom:
            raise ValueError("Chromosome cannot be empty.")
        if not self.side:
            raise ValueError("Side cannot be empty.")
        if isinstance(self.quirk_type, str) and not isinstance(self.quirk_type, IntervalQuirk):
            self.quirk_type = IntervalQuirk(self.quirk_type)
        if not isinstance(self.unitig_ids, tuple):
            self.unitig_ids = tuple(self.unitig_ids)
        if not isinstance(self.excluded_unitig_ids, tuple):
            self.excluded_unitig_ids = tuple(self.excluded_unitig_ids)
        if (
            self.length_min_bp is not None
            and self.length_max_bp is not None
            and self.length_min_bp > self.length_max_bp
        ):
            raise ValueError(
                f"length_min_bp ({self.length_min_bp}) cannot exceed length_max_bp ({self.length_max_bp})."
            )
        if self.agp_start is not None and self.agp_end is not None and self.agp_start > self.agp_end:
            raise ValueError(f"agp_start ({self.agp_start}) cannot exceed agp_end ({self.agp_end}).")

    def to_dict(self) -> dict[str, Any]:
        """Serialize record to dictionary with string enum representations."""
        d = asdict(self)
        d["quirk_type"] = self.quirk_type.value
        d["unitig_ids"] = list(self.unitig_ids)
        d["excluded_unitig_ids"] = list(self.excluded_unitig_ids)
        return d

    def to_tsv_row(self) -> list[str]:
        """Serialize record to an ordered list of string values for TSV emission."""
        return [
            self.feature_id,
            self.chrom,
            self.side,
            self.quirk_type.value,
            str(self.spanning_read_count) if self.spanning_read_count is not None else "",
            str(self.length_min_bp) if self.length_min_bp is not None else "",
            str(self.length_median_bp) if self.length_median_bp is not None else "",
            str(self.length_max_bp) if self.length_max_bp is not None else "",
            str(self.length_iqr_bp) if self.length_iqr_bp is not None else "",
            str(self.repeat_monomer_len_bp) if self.repeat_monomer_len_bp is not None else "",
            self.copy_number_range or "",
            self.left_anchor_port or "",
            self.right_anchor_port or "",
            self.agp_component_type or "",
            self.agp_linkage_evidence or "",
            ",".join(self.unitig_ids),
            str(self.read_count) if self.read_count is not None else "",
            str(self.core_coverage_bp_median) if self.core_coverage_bp_median is not None else "",
            str(self.core_depth_ratio) if self.core_depth_ratio is not None else "",
            str(self.representative_sequence_selected),
            self.sequence_provenance or "",
            str(self.unpolished),
            self.notes,
            str(self.agp_start) if self.agp_start is not None else "",
            str(self.agp_end) if self.agp_end is not None else "",
            ",".join(self.excluded_unitig_ids),
        ]


QUIRK_TSV_COLUMNS: list[str] = [
    "feature_id",
    "chrom",
    "side",
    "quirk_type",
    "spanning_read_count",
    "length_min_bp",
    "length_median_bp",
    "length_max_bp",
    "length_iqr_bp",
    "repeat_monomer_len_bp",
    "copy_number_range",
    "left_anchor_port",
    "right_anchor_port",
    "agp_component_type",
    "agp_linkage_evidence",
    "unitig_ids",
    "read_count",
    "core_coverage_bp_median",
    "core_depth_ratio",
    "representative_sequence_selected",
    "sequence_provenance",
    "unpolished",
    "notes",
    "agp_start",
    "agp_end",
    "excluded_unitig_ids",
]

QUIRK_COLUMNS: list[str] = QUIRK_TSV_COLUMNS
RESOLUTION_COLUMNS: list[str] = QUIRK_COLUMNS
RESOLUTION_TSV_COLUMNS: list[str] = QUIRK_TSV_COLUMNS

ResolutionType = IntervalQuirk
ResolutionRecord = QuirkRecord


def _clean_tsv_stream(stream: Iterable[str]) -> list[str]:
    """Filter comment and blank lines, normalizing commented headers."""
    cleaned: list[str] = []
    for line in stream:
        s = line.strip()
        if not s:
            continue
        if s.startswith("#") and "\t" in s and s.lstrip("#").strip().startswith("feature_id"):
            cleaned.append(s.lstrip("#").strip() + "\n")
            continue
        if s.startswith("#"):
            continue
        cleaned.append(line)
    return cleaned


def emit_quirks_tsv(
    records: Iterable[QuirkRecord],
    out_file: TextIO | Path | str | None = None,
) -> str:
    """Emit a list or iterable of QuirkRecords as a standardized TSV string (and optional file).

    Args:
        records: Iterable of QuirkRecord objects.
        out_file: Optional destination file path or open text file object.

    Returns:
        Formatted TSV content as a string.
    """
    buffer = StringIO()
    writer = csv.writer(buffer, delimiter="\t", lineterminator="\n")
    writer.writerow(QUIRK_TSV_COLUMNS)
    for rec in records:
        writer.writerow(rec.to_tsv_row())
    content = buffer.getvalue()

    if out_file is not None:
        if isinstance(out_file, str | Path):
            with open(out_file, "w", encoding="utf-8") as f:
                f.write(content)
        else:
            out_file.write(content)

    return content


emit_resolutions_tsv = emit_quirks_tsv


def _to_int_or_none(row: dict[str, str], key: str) -> int | None:
    val = row.get(key, "").strip()
    return int(val) if val else None


def _to_float_or_none(row: dict[str, str], key: str) -> float | None:
    val = row.get(key, "").strip()
    return float(val) if val else None


def parse_quirks_tsv(file: TextIO | Path | str) -> list[QuirkRecord]:
    """Parse a standardized quirks or resolutions TSV content string, stream, or file into a list of QuirkRecord objects.

    Args:
        file: Input file path, open text stream, or TSV content string.

    Returns:
        List of deserialized QuirkRecord instances.
    """

    def _parse_rows(reader: csv.DictReader[str]) -> list[QuirkRecord]:
        records: list[QuirkRecord] = []
        for row in reader:
            if not row:
                continue
            fid = row.get("feature_id")
            if not fid or fid.startswith("#"):
                continue
            utg_str = row.get("unitig_ids", "")
            unitig_ids = tuple(x.strip() for x in utg_str.split(",") if x.strip())
            excl_utg_str = row.get("excluded_unitig_ids", "")
            excluded_unitig_ids = tuple(x.strip() for x in excl_utg_str.split(",") if x.strip())
            rep_sel_raw = row.get("representative_sequence_selected", "").lower()
            rep_sel = rep_sel_raw in ("true", "1")
            seq_prov = row.get("sequence_provenance", "").strip() or None
            unpol_raw = row.get("unpolished", "").strip().lower()
            unpolished = unpol_raw in ("true", "1")

            records.append(
                QuirkRecord(
                    feature_id=fid.strip(),
                    chrom=row["chrom"].strip(),
                    side=row["side"].strip(),
                    quirk_type=IntervalQuirk(row["quirk_type"].strip()),
                    spanning_read_count=_to_int_or_none(row, "spanning_read_count"),
                    length_min_bp=_to_int_or_none(row, "length_min_bp"),
                    length_median_bp=_to_int_or_none(row, "length_median_bp"),
                    length_max_bp=_to_int_or_none(row, "length_max_bp"),
                    length_iqr_bp=_to_int_or_none(row, "length_iqr_bp"),
                    repeat_monomer_len_bp=_to_int_or_none(row, "repeat_monomer_len_bp"),
                    copy_number_range=row.get("copy_number_range") or None,
                    left_anchor_port=row.get("left_anchor_port") or None,
                    right_anchor_port=row.get("right_anchor_port") or None,
                    agp_component_type=row.get("agp_component_type") or None,
                    agp_linkage_evidence=row.get("agp_linkage_evidence") or None,
                    unitig_ids=unitig_ids,
                    read_count=_to_int_or_none(row, "read_count"),
                    core_coverage_bp_median=_to_float_or_none(row, "core_coverage_bp_median"),
                    core_depth_ratio=_to_float_or_none(row, "core_depth_ratio"),
                    representative_sequence_selected=rep_sel,
                    sequence_provenance=seq_prov,
                    unpolished=unpolished,
                    notes=row.get("notes", ""),
                    agp_start=_to_int_or_none(row, "agp_start"),
                    agp_end=_to_int_or_none(row, "agp_end"),
                    excluded_unitig_ids=excluded_unitig_ids,
                )
            )
        return records

    if isinstance(file, Path):
        with open(file, encoding="utf-8") as f:
            return _parse_rows(csv.DictReader(_clean_tsv_stream(f), delimiter="\t"))
    if isinstance(file, str):
        if "\n" not in file and "\t" not in file:
            stripped = file.strip()
            if not stripped:
                return []
            p = Path(file)
            if p.suffix in (".tsv", ".txt", ".tab") or "/" in file or "\\" in file or p.is_file():
                with open(p, encoding="utf-8") as f:
                    return _parse_rows(csv.DictReader(_clean_tsv_stream(f), delimiter="\t"))
            return _parse_rows(csv.DictReader(StringIO(file), delimiter="\t"))
        return _parse_rows(csv.DictReader(_clean_tsv_stream(StringIO(file)), delimiter="\t"))
    return _parse_rows(csv.DictReader(_clean_tsv_stream(file), delimiter="\t"))


parse_resolutions_tsv = parse_quirks_tsv


END_STATE_TSV_COLUMNS: list[str] = [
    "feature_id",
    "chrom",
    "side",
    "state",
    "somatic_port",
    "germline_port",
    "shared_port",
    "core_coverage_bp_median",
    "somatic_depth_ratio",
    "germline_depth_ratio",
    "resolved",
    "notes",
]

JUNCTION_TSV_COLUMNS: list[str] = [
    "feature_id",
    "chrom",
    "side",
    "junction_class",
    "stem_port",
    "somatic_terminus_port",
    "germline_continuation_port",
    "break_index",
    "core_depth_ratio",
    "read_count",
    "confidence_tier",
    "notes",
]


def emit_end_states_tsv(
    records: Sequence[EndStateRecord],
    out_file: TextIO | Path | str | None = None,
) -> str:
    """Emit a sequence of EndStateRecords as a standardized TSV string (and optional file).

    Args:
        records: Sequence of EndStateRecord objects.
        out_file: Optional destination file path or open text file object.

    Returns:
        Formatted TSV content as a string.
    """
    buffer = StringIO()
    writer = csv.writer(buffer, delimiter="\t", lineterminator="\n")
    writer.writerow(END_STATE_TSV_COLUMNS)
    for rec in records:
        writer.writerow(rec.to_tsv_row())
    content = buffer.getvalue()

    if out_file is not None:
        if isinstance(out_file, str | Path):
            with open(out_file, "w", encoding="utf-8") as f:
                f.write(content)
        else:
            out_file.write(content)

    return content


def parse_end_states_tsv(tsv_content: str | Path | TextIO) -> list[EndStateRecord]:
    """Parse a standardized end_states TSV content string or file into a list of EndStateRecords.

    Args:
        tsv_content: TSV content string, open text stream, or file path.

    Returns:
        List of deserialized EndStateRecord instances.
    """

    def _parse_reader(reader: csv.DictReader[str]) -> list[EndStateRecord]:
        records: list[EndStateRecord] = []
        for row in reader:
            if not row:
                continue
            fid = row.get("feature_id")
            if not fid or fid.startswith("#"):
                continue
            resolved_raw = row.get("resolved", "").strip().lower()
            resolved = resolved_raw in ("true", "1", "yes", "t") if resolved_raw else True

            records.append(
                EndStateRecord(
                    feature_id=fid.strip(),
                    chrom=row["chrom"].strip(),
                    side=row["side"].strip(),
                    state=EndState(row["state"].strip()),
                    somatic_port=(row.get("somatic_port") or "").strip() or None,
                    germline_port=(row.get("germline_port") or "").strip() or None,
                    shared_port=(row.get("shared_port") or "").strip() or None,
                    core_coverage_bp_median=_to_float_or_none(row, "core_coverage_bp_median"),
                    somatic_depth_ratio=_to_float_or_none(row, "somatic_depth_ratio"),
                    germline_depth_ratio=_to_float_or_none(row, "germline_depth_ratio"),
                    resolved=resolved,
                    notes=row.get("notes", ""),
                )
            )
        return records

    if isinstance(tsv_content, Path):
        with open(tsv_content, encoding="utf-8") as f:
            return _parse_reader(csv.DictReader(_clean_tsv_stream(f), delimiter="\t"))
    if isinstance(tsv_content, str):
        if "\n" not in tsv_content and "\t" not in tsv_content:
            stripped = tsv_content.strip()
            if not stripped:
                return []
            p = Path(tsv_content)
            if p.suffix in (".tsv", ".txt", ".tab") or "/" in tsv_content or "\\" in tsv_content or p.is_file():
                with open(p, encoding="utf-8") as f:
                    return _parse_reader(csv.DictReader(_clean_tsv_stream(f), delimiter="\t"))
            return _parse_reader(csv.DictReader(StringIO(tsv_content), delimiter="\t"))
        return _parse_reader(csv.DictReader(_clean_tsv_stream(StringIO(tsv_content)), delimiter="\t"))
    return _parse_reader(csv.DictReader(_clean_tsv_stream(tsv_content), delimiter="\t"))


def emit_junctions_tsv(
    records: Sequence[JunctionRecord],
    out_file: TextIO | Path | str | None = None,
) -> str:
    """Emit a sequence of JunctionRecords as a standardized TSV string (and optional file).

    Args:
        records: Sequence of JunctionRecord objects.
        out_file: Optional destination file path or open text file object.

    Returns:
        Formatted TSV content as a string.
    """
    buffer = StringIO()
    writer = csv.writer(buffer, delimiter="\t", lineterminator="\n")
    writer.writerow(JUNCTION_TSV_COLUMNS)
    for rec in records:
        writer.writerow(rec.to_tsv_row())
    content = buffer.getvalue()

    if out_file is not None:
        if isinstance(out_file, str | Path):
            with open(out_file, "w", encoding="utf-8") as f:
                f.write(content)
        else:
            out_file.write(content)

    return content


def parse_junctions_tsv(tsv_content: str | Path | TextIO) -> list[JunctionRecord]:
    """Parse a standardized junctions TSV content string or file into a list of JunctionRecords.

    Args:
        tsv_content: TSV content string, open text stream, or file path.

    Returns:
        List of deserialized JunctionRecord instances.
    """

    def _parse_reader(reader: csv.DictReader[str]) -> list[JunctionRecord]:
        records: list[JunctionRecord] = []
        for row in reader:
            if not row:
                continue
            fid = row.get("feature_id")
            if not fid or fid.startswith("#"):
                continue
            records.append(
                JunctionRecord(
                    feature_id=fid.strip(),
                    chrom=row["chrom"].strip(),
                    side=row["side"].strip(),
                    junction_class=JunctionClass(row["junction_class"].strip()),
                    stem_port=row["stem_port"].strip(),
                    somatic_terminus_port=(row.get("somatic_terminus_port") or "").strip() or None,
                    germline_continuation_port=(row.get("germline_continuation_port") or "").strip() or None,
                    break_index=int(row.get("break_index", 1) or 1),
                    core_depth_ratio=_to_float_or_none(row, "core_depth_ratio"),
                    read_count=_to_int_or_none(row, "read_count"),
                    confidence_tier=row.get("confidence_tier") or "CONFIDENT",
                    notes=row.get("notes", ""),
                )
            )
        return records

    if isinstance(tsv_content, Path):
        with open(tsv_content, encoding="utf-8") as f:
            return _parse_reader(csv.DictReader(_clean_tsv_stream(f), delimiter="\t"))
    if isinstance(tsv_content, str):
        if "\n" not in tsv_content and "\t" not in tsv_content:
            stripped = tsv_content.strip()
            if not stripped:
                return []
            p = Path(tsv_content)
            if p.suffix in (".tsv", ".txt", ".tab") or "/" in tsv_content or "\\" in tsv_content or p.is_file():
                with open(p, encoding="utf-8") as f:
                    return _parse_reader(csv.DictReader(_clean_tsv_stream(f), delimiter="\t"))
            return _parse_reader(csv.DictReader(StringIO(tsv_content), delimiter="\t"))
        return _parse_reader(csv.DictReader(_clean_tsv_stream(StringIO(tsv_content)), delimiter="\t"))
    return _parse_reader(csv.DictReader(_clean_tsv_stream(tsv_content), delimiter="\t"))


def create_junction_record(
    call: PdeJunctionCall,
    chrom: str,
    side: str,
    junction_class: JunctionClass | str,
    confidence_tier: str = "CONFIDENT",
    break_index: int = 1,
    read_count: int | None = None,
    notes: str = "",
) -> JunctionRecord:
    """Construct a Layer 2 JunctionRecord from a classifier PdeJunctionCall.

    Args:
        call: PdeJunctionCall output from PDE junction classification.
        chrom: Chromosome name (e.g. 'I', 'IV', 'X').
        side: Arm end or locus ('5p', '3p', 'mid').
        junction_class: JunctionClass enum or value ('pde_terminal' or 'pde_intercalary').
        confidence_tier: Classification confidence tier (default 'CONFIDENT').
        break_index: 1-based index for intercalary breaks.
        read_count: Optional supporting read count across the junction.
        notes: Contextual notes or provenance information.

    Returns:
        Canonical JunctionRecord instance.
    """
    j_class = JunctionClass(junction_class) if not isinstance(junction_class, JunctionClass) else junction_class

    if j_class == JunctionClass.PDE_INTERCALARY:
        fid = build_feature_id(j_class, chrom, side, index=break_index)
    else:
        fid = build_feature_id(j_class, chrom, side, index=break_index if break_index > 1 else None)

    stem_port = call.port.to_symbol() if hasattr(call.port, "to_symbol") else str(call.port)
    somatic_port = (
        call.somatic_port.to_symbol()
        if hasattr(call.somatic_port, "to_symbol")
        else (str(call.somatic_port) if call.somatic_port is not None else None)
    )
    germline_port = (
        call.germline_port.to_symbol()
        if hasattr(call.germline_port, "to_symbol")
        else (str(call.germline_port) if call.germline_port is not None else None)
    )

    final_notes = notes
    if getattr(call, "is_repeat_entry", False):
        final_notes = f"{notes}; repeat_entry" if notes else "repeat_entry"

    core_depth_ratio: float | None = None
    if getattr(call, "core_depth_ratio", 0.0) != 0.0:
        core_depth_ratio = call.core_depth_ratio
    elif getattr(call, "depth_ratio", 0.0) != 0.0:
        core_depth_ratio = call.depth_ratio
    elif getattr(call, "core_depth_ratio", None) is not None:
        core_depth_ratio = call.core_depth_ratio
    elif getattr(call, "depth_ratio", None) is not None:
        core_depth_ratio = call.depth_ratio

    return JunctionRecord(
        feature_id=fid,
        chrom=chrom,
        side=side,
        junction_class=j_class,
        stem_port=stem_port,
        somatic_terminus_port=somatic_port,
        germline_continuation_port=germline_port,
        break_index=break_index,
        core_depth_ratio=core_depth_ratio,
        read_count=read_count,
        confidence_tier=confidence_tier,
        notes=final_notes,
    )


def _reaches_germline_telomere_cap(
    graph: BidirectedGraph,
    ledger_by_port: Mapping[Any, EndEvidence],
    start_port: SegmentPort,
    max_hops: int,
) -> list[SegmentPort] | None:
    """Bounded walk from a germline entry port to identify an authentic terminal telomeric cap."""

    def _is_contaminant(name: str) -> bool:
        contam_classes = ("bacteria", "archaea", "prokarya", "contaminant")
        for key in (f"{name}+", f"{name}-", f"{name}:5p", f"{name}:3p", name):
            ev = ledger_by_port.get(key)
            if ev is not None:
                if getattr(ev, "is_contaminant", False):
                    return True
                cls_str = str(getattr(ev, "tiara_class", "")).lower().strip()
                if cls_str in contam_classes:
                    return True
        return False

    def _get_usco(name: str) -> int:
        for key in (f"{name}+", f"{name}-", f"{name}:5p", f"{name}:3p", name):
            ev = ledger_by_port.get(key)
            if ev is not None and getattr(ev, "usco_hit_count", 0) > 0:
                return ev.usco_hit_count
        return 0

    ev_start = _lookup_evidence(ledger_by_port, start_port)
    if ev_start is not None and ev_start.has_terminal_telomere:
        deg = graph.degree(start_port) if start_port.name in graph.segments else ev_start.degree
        if deg == 0:
            return [start_port]

    if start_port.name not in graph.segments:
        return None

    queue: list[tuple[SegmentPort, list[SegmentPort]]] = [(start_port, [start_port])]
    visited: set[SegmentPort] = {start_port}

    while queue:
        curr, path = queue.pop(0)
        if len(path) - 1 >= max_hops:
            continue

        neighbors = (
            list(graph.neighbors(curr))
            if hasattr(graph, "neighbors")
            else [e.to_port for e in graph.get_outgoing(curr)]
        )

        for nxt in neighbors:
            if nxt in visited:
                continue
            if _get_usco(nxt.name) != 0 or _is_contaminant(nxt.name):
                continue

            nxt_ev = _lookup_evidence(ledger_by_port, nxt)
            nxt_deg = graph.degree(nxt) if nxt.name in graph.segments else (nxt_ev.degree if nxt_ev else 0)
            if nxt_ev is not None and nxt_ev.has_terminal_telomere and nxt_deg == 0:
                return [*path, nxt]

            visited.add(nxt)
            queue.append((nxt, [*path, nxt]))

    return None


def collect_features(
    graph: BidirectedGraph,
    ledger: Sequence[EndEvidence] | Mapping[Any, EndEvidence],
    core_coverage: float,
    policy: Policy = DEFAULT_POLICY,
) -> tuple[list[JunctionRecord], list[EndStateRecord], list[QuirkRecord]]:
    """Aggregate Layer 1 EndStates, Layer 2 Junctions, and Layer 3 Quirks from graph evidence."""
    ledger_list: list[EndEvidence]
    ledger_by_port: dict[Any, EndEvidence] = {}
    if isinstance(ledger, Sequence) and not isinstance(ledger, Mapping):
        ledger_list = list(ledger)
        for ev in ledger_list:
            ledger_by_port[ev.port] = ev
            ledger_by_port[ev.port.to_symbol()] = ev
            ledger_by_port[f"{ev.segment_name}:{ev.end}"] = ev
    else:
        ledger_by_port = dict(ledger)
        ledger_list = [ev for ev in ledger.values() if isinstance(ev, EndEvidence)]

    ledger_list.sort(key=lambda e: (e.segment_name, e.end))

    junction_records: list[JunctionRecord] = []
    end_state_records: list[EndStateRecord] = []
    quirk_records: list[QuirkRecord] = []

    accounted_ports: set[str] = set()

    # Step 1: Detect PDE junctions across all candidate ports
    pde_calls_with_ev: list[tuple[EndEvidence, PdeJunctionCall]] = []
    for ev in ledger_list:
        pde_call = classify_pde_junction(
            ev.port,
            graph,
            ledger_by_port,
            policy=policy,
            core_coverage=core_coverage,
        )
        if pde_call.is_pde_junction:
            pde_calls_with_ev.append((ev, pde_call))

    called_stem_ports: set[Any] = set()
    for _, c in pde_calls_with_ev:
        sp = getattr(c, "stem_port", getattr(c, "port", None))
        if sp is not None:
            called_stem_ports.add(sp)
            if hasattr(sp, "to_symbol"):
                called_stem_ports.add(sp.to_symbol())
            if hasattr(sp, "name"):
                called_stem_ports.add(sp.name)
            if isinstance(sp, str):
                called_stem_ports.add(graph._resolve_port(sp))

    def _get_seg_usco(s_name: str) -> int:
        for k in (f"{s_name}+", f"{s_name}-", f"{s_name}:5p", f"{s_name}:3p", s_name):
            ev_item = ledger_by_port.get(k)
            if ev_item is not None and getattr(ev_item, "usco_hit_count", 0) > 0:
                return ev_item.usco_hit_count
        return 0

    for ev, pde_call in pde_calls_with_ev:
        j_class = JunctionClass.PDE_TERMINAL if pde_call.somatic_port is not None else JunctionClass.PDE_INTERCALARY
        j_rec = create_junction_record(
            call=pde_call,
            chrom=ev.segment_name,
            side=ev.end,
            junction_class=j_class,
        )
        junction_records.append(j_rec)

        accounted_ports.add(ev.port.to_symbol())
        if pde_call.somatic_port:
            accounted_ports.add(pde_call.somatic_port.to_symbol())
        if pde_call.germline_port:
            accounted_ports.add(pde_call.germline_port.to_symbol())

        is_som_dropout = getattr(pde_call, "is_somatic_dropout", False)
        som_capped = (not is_som_dropout) and (pde_call.somatic_port is not None)
        germ_cap_path: list[SegmentPort] | None = None
        if pde_call.germline_port is not None:
            max_bubble_hops = (
                policy.graph.max_bubble_search_hops.value
                if hasattr(policy, "graph") and hasattr(policy.graph, "max_bubble_search_hops")
                else policy.pde.max_somatic_cap_hops.value
            )
            germ_cap_path = _reaches_germline_telomere_cap(
                graph=graph,
                ledger_by_port=ledger_by_port,
                start_port=pde_call.germline_port,
                max_hops=max_bubble_hops,
            )
        germ_capped = germ_cap_path is not None
        if germ_cap_path:
            for gp in germ_cap_path:
                accounted_ports.add(gp.to_symbol())

        germ_intercalary = False
        if pde_call.germline_port is not None and pde_call.germline_port.name in graph.segments and not germ_capped:
            if hasattr(graph, "neighbors"):
                neighbors = list(graph.neighbors(pde_call.germline_port))
            else:
                neighbors = [e.to_port for e in graph.get_outgoing(pde_call.germline_port)]

            for nb in neighbors:
                inv_nb = nb.invert() if hasattr(nb, "invert") else (nb.inverted() if hasattr(nb, "inverted") else ~nb)
                if (
                    inv_nb in called_stem_ports
                    or nb in called_stem_ports
                    or (hasattr(inv_nb, "to_symbol") and inv_nb.to_symbol() in called_stem_ports)
                    or (hasattr(nb, "to_symbol") and nb.to_symbol() in called_stem_ports)
                ):
                    germ_intercalary = True
                    break
                nb_name = nb.name if hasattr(nb, "name") else str(nb).rstrip("+-")
                if nb_name != ev.segment_name and _get_seg_usco(nb_name) > 0:
                    germ_intercalary = True
                    break

        if som_capped and (germ_capped or germ_intercalary):
            state = EndState.END_COMPLETE
            resolved = True
        elif som_capped and not (germ_capped or germ_intercalary):
            state = EndState.END_GERMLINE_ORPHAN
            resolved = False
        elif not som_capped and (germ_capped or germ_intercalary):
            state = EndState.END_SOMATIC_MISSING
            resolved = False
        else:
            state = EndState.END_SOMATIC_MISSING if is_som_dropout else EndState.END_GERMLINE_ORPHAN
            resolved = False

        fid = build_feature_id(state, ev.segment_name, ev.end)
        som_sym = pde_call.somatic_port.to_symbol() if pde_call.somatic_port else None
        germ_sym = pde_call.germline_port.to_symbol() if pde_call.germline_port else None

        som_seg = graph.segments.get(pde_call.somatic_port.name) if pde_call.somatic_port else None
        som_depth = som_seg.depth if som_seg is not None else None
        som_ratio = (som_depth / core_coverage) if (som_depth is not None and core_coverage > 0) else None
        germ_ratio = (
            pde_call.core_depth_ratio
            if pde_call.core_depth_ratio > 0
            else (pde_call.depth_ratio if pde_call.depth_ratio > 0 else None)
        )

        end_rec = EndStateRecord(
            feature_id=fid,
            chrom=ev.segment_name,
            side=ev.end,
            state=state,
            somatic_port=som_sym,
            germline_port=germ_sym,
            core_coverage_bp_median=core_coverage,
            somatic_depth_ratio=som_ratio,
            germline_depth_ratio=germ_ratio,
            resolved=resolved,
            notes=f"PDE junction at {ev.port.to_symbol()}",
        )
        end_state_records.append(end_rec)

        if not (germ_capped or germ_intercalary):
            q_fid = build_feature_id(IntervalQuirk.SAT_GAP_VARIABLE, ev.segment_name, ev.end, index=1)
            q_rec = QuirkRecord(
                feature_id=q_fid,
                chrom=ev.segment_name,
                side=ev.end,
                quirk_type=IntervalQuirk.SAT_GAP_VARIABLE,
                left_anchor_port=germ_sym or ev.port.to_symbol(),
                agp_component_type="U",
                agp_linkage_evidence="unbridged",
                unitig_ids=(ev.segment_name,),
                core_coverage_bp_median=core_coverage,
                representative_sequence_selected=False,
                notes="Unbridged variable satellite gap",
            )
            quirk_records.append(q_rec)

    # Step 2: Evaluate remaining ports for capped ends, orphan caps, and tips
    for ev in ledger_list:
        sym = ev.port.to_symbol()
        if sym in accounted_ports:
            continue

        deg = graph.degree(ev.port) if ev.port.name in graph.segments else ev.degree
        if deg == 0 and ev.has_terminal_telomere:
            other_strand = Strand.FORWARD if ev.port.strand == Strand.MINUS else Strand.MINUS
            other_port = SegmentPort(ev.segment_name, other_strand)
            other_deg = graph.degree(other_port) if ev.segment_name in graph.segments else 0
            other_ev = _lookup_evidence(ledger_by_port, other_port)

            if other_deg == 0 and (other_ev is None or not other_ev.has_terminal_telomere):
                # Orphan cap: inward end is open, telomere end is capped
                pass
            else:
                if ev.rdna_hit_count > 0 or ev.has_terminal_rdna:
                    state = EndState.END_SHARED_EXEMPT
                    q_fid = build_feature_id(IntervalQuirk.RDNA_CASSETTE, ev.segment_name, ev.end, index=1)
                    quirk_records.append(
                        QuirkRecord(
                            feature_id=q_fid,
                            chrom=ev.segment_name,
                            side=ev.end,
                            quirk_type=IntervalQuirk.RDNA_CASSETTE,
                            left_anchor_port=ev.port.to_symbol(),
                            agp_component_type="W",
                            agp_linkage_evidence="graph",
                            unitig_ids=(ev.segment_name,),
                            notes="rDNA cluster at chromosome terminus",
                        )
                    )
                else:
                    state = EndState.END_COMPLETE

                fid = build_feature_id(state, ev.segment_name, ev.end)
                end_state_records.append(
                    EndStateRecord(
                        feature_id=fid,
                        chrom=ev.segment_name,
                        side=ev.end,
                        state=state,
                        shared_port=ev.port.to_symbol(),
                        core_coverage_bp_median=core_coverage,
                        resolved=True,
                        notes="Terminal telomeric cap",
                    )
                )
        elif deg == 0 and not ev.has_terminal_telomere:
            state = EndState.END_GERMLINE_ORPHAN
            fid = build_feature_id(state, ev.segment_name, ev.end)
            end_state_records.append(
                EndStateRecord(
                    feature_id=fid,
                    chrom=ev.segment_name,
                    side=ev.end,
                    state=state,
                    germline_port=ev.port.to_symbol(),
                    core_coverage_bp_median=core_coverage,
                    resolved=False,
                    notes="Uncapped terminal tip",
                )
            )

    junction_records.sort(key=lambda r: r.feature_id)
    end_state_records.sort(key=lambda r: r.feature_id)
    quirk_records.sort(key=lambda r: r.feature_id)
    return junction_records, end_state_records, quirk_records
