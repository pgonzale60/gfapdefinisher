"""Biological quirk, junction, and chromosome end-state taxonomy.

Defines the ratified 3-layer biological feature schema for assembly graph
finishing, capturing chromosome end completeness, somatic diminution (PDE)
junctions, and interval-level structural quirks.

Ratified in Issue #26.
"""

from __future__ import annotations

import csv
from collections.abc import Iterable
from dataclasses import asdict, dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any, TextIO

__all__ = [
    "EndState",
    "JunctionClass",
    "IntervalQuirk",
    "EndStateRecord",
    "JunctionRecord",
    "QuirkRecord",
    "build_feature_id",
    "emit_quirks_tsv",
    "parse_quirks_tsv",
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

    # Lowercase aliases
    sat_gap_variable = "sat_gap_variable"
    sat_gap_fixed = "sat_gap_fixed"
    rdna_cassette = "rdna_cassette"
    uncollapsed_copy = "uncollapsed_copy"
    low_support_stub = "low_support_stub"
    bubble_rejoining = "bubble_rejoining"
    itr = "itr"


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
        core_depth_ratio: Depth ratio of germline branch relative to somatic stem.
        read_count: Supporting read count across the junction boundary.
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
    notes: str = ""

    def __post_init__(self) -> None:
        if not self.chrom:
            raise ValueError("Chromosome cannot be empty.")
        if not self.side:
            raise ValueError("Side cannot be empty.")
        if not self.stem_port:
            raise ValueError("Stem port cannot be empty.")
        if isinstance(self.junction_class, str) and not isinstance(self.junction_class, JunctionClass):
            self.junction_class = JunctionClass(self.junction_class)

    def to_dict(self) -> dict[str, Any]:
        """Serialize record to dictionary with string enum representations."""
        d = asdict(self)
        d["junction_class"] = self.junction_class.value
        return d


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
    notes: str = ""

    def __post_init__(self) -> None:
        if not self.chrom:
            raise ValueError("Chromosome cannot be empty.")
        if not self.side:
            raise ValueError("Side cannot be empty.")
        if isinstance(self.quirk_type, str) and not isinstance(self.quirk_type, IntervalQuirk):
            self.quirk_type = IntervalQuirk(self.quirk_type)
        if not isinstance(self.unitig_ids, tuple):
            self.unitig_ids = tuple(self.unitig_ids)
        if (
            self.length_min_bp is not None
            and self.length_max_bp is not None
            and self.length_min_bp > self.length_max_bp
        ):
            raise ValueError(
                f"length_min_bp ({self.length_min_bp}) cannot exceed length_max_bp ({self.length_max_bp})."
            )

    def to_dict(self) -> dict[str, Any]:
        """Serialize record to dictionary with string enum representations."""
        d = asdict(self)
        d["quirk_type"] = self.quirk_type.value
        d["unitig_ids"] = list(self.unitig_ids)
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
            self.notes,
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
    "notes",
]


def emit_quirks_tsv(records: Iterable[QuirkRecord], out_file: TextIO | Path | str) -> None:
    """Emit a list or iterable of QuirkRecords as a standardized TSV file.

    Args:
        records: Iterable of QuirkRecord objects.
        out_file: Destination file path or open text file object.
    """

    def _write_records(writer: Any) -> None:
        writer.write("\t".join(QUIRK_TSV_COLUMNS) + "\n")
        for rec in records:
            writer.write("\t".join(rec.to_tsv_row()) + "\n")

    if isinstance(out_file, str | Path):
        with open(out_file, "w", encoding="utf-8") as f:
            _write_records(f)
    else:
        _write_records(out_file)


def _to_int_or_none(row: dict[str, str], key: str) -> int | None:
    val = row.get(key, "").strip()
    return int(val) if val else None


def _to_float_or_none(row: dict[str, str], key: str) -> float | None:
    val = row.get(key, "").strip()
    return float(val) if val else None


def parse_quirks_tsv(file: TextIO | Path | str) -> list[QuirkRecord]:
    """Parse a standardized quirks.tsv file into a list of QuirkRecord objects.

    Args:
        file: Input file path or open text stream.

    Returns:
        List of deserialized QuirkRecord instances.
    """

    def _parse_rows(reader: csv.DictReader[str]) -> list[QuirkRecord]:
        records: list[QuirkRecord] = []
        for row in reader:
            utg_str = row.get("unitig_ids", "")
            unitig_ids = tuple(x.strip() for x in utg_str.split(",") if x.strip())
            rep_sel_raw = row.get("representative_sequence_selected", "").lower()
            rep_sel = rep_sel_raw in ("true", "1")

            records.append(
                QuirkRecord(
                    feature_id=row["feature_id"],
                    chrom=row["chrom"],
                    side=row["side"],
                    quirk_type=IntervalQuirk(row["quirk_type"]),
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
                    notes=row.get("notes", ""),
                )
            )
        return records

    if isinstance(file, str | Path):
        with open(file, encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            return _parse_rows(reader)
    reader = csv.DictReader(file, delimiter="\t")
    return _parse_rows(reader)
