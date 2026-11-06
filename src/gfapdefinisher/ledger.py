"""Unitig end-evidence ledger generator for assembly graph nodes.

Provides data structures and functions to extract Tier B evidence for both
5' and 3' ports of every unitig in a bidirected assembly graph:
- Graph degree, tip status, and outgoing neighbor ports
- TIDK telomeric repeat counts within the terminal window
- Barrnap ribosomal RNA hits within the terminal window
- Compleasm universal single-copy ortholog (USCO) gene hits
- Tiara taxonomic/contaminant classifications
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, TextIO

from gfapdefinisher.annotations import NodeAnnotations
from gfapdefinisher.graph import BidirectedGraph
from gfapdefinisher.models import SegmentPort, Strand
from gfapdefinisher.policy import DEFAULT_POLICY, Policy
from gfapdefinisher.telomere import scan_terminal_tract

DEFAULT_TERMINAL_WINDOW_BP: int = DEFAULT_POLICY.telomere.terminal_scan_bp.value


@dataclass(slots=True, frozen=True)
class EndEvidence:
    """Evidence record for one oriented end of an assembly segment.

    Captures graph topology, outgoing dovetails, repeat evidence,
    and gene content within the terminal window.
    """

    segment_name: str
    end: str
    port: SegmentPort
    length: int
    depth: float
    degree: int
    is_tip: bool
    outgoing_neighbors: tuple[str, ...]
    telomere_canonical_count: int
    telomere_opposite_count: int
    has_terminal_telomere: bool
    rdna_hit_count: int
    rdna_features: tuple[str, ...]
    has_terminal_rdna: bool
    usco_hit_count: int
    usco_genes: tuple[str, ...]
    tiara_class: str
    telomere_fwd_count: int = 0
    telomere_rev_count: int = 0
    telomere_offset_bp: int | None = None
    telomere_tract_bp: int | None = None

    def __init__(
        self,
        segment_name: str,
        end: str,
        port: SegmentPort,
        length: int,
        depth: float,
        degree: int,
        is_tip: bool,
        outgoing_neighbors: tuple[str, ...],
        telomere_canonical_count: int | None = None,
        telomere_opposite_count: int | None = None,
        has_terminal_telomere: bool = False,
        rdna_hit_count: int = 0,
        rdna_features: tuple[str, ...] = (),
        has_terminal_rdna: bool = False,
        usco_hit_count: int = 0,
        usco_genes: tuple[str, ...] = (),
        tiara_class: str = "",
        telomere_fwd_count: int | None = None,
        telomere_rev_count: int | None = None,
        telomere_offset_bp: int | None = None,
        telomere_tract_bp: int | None = None,
    ) -> None:
        object.__setattr__(self, "segment_name", segment_name)
        object.__setattr__(self, "end", end)
        object.__setattr__(self, "port", port)
        object.__setattr__(self, "length", length)
        object.__setattr__(self, "depth", depth)
        object.__setattr__(self, "degree", degree)
        object.__setattr__(self, "is_tip", is_tip)
        object.__setattr__(self, "outgoing_neighbors", outgoing_neighbors)
        object.__setattr__(self, "has_terminal_telomere", has_terminal_telomere)
        object.__setattr__(self, "rdna_hit_count", rdna_hit_count)
        object.__setattr__(self, "rdna_features", rdna_features)
        object.__setattr__(self, "has_terminal_rdna", has_terminal_rdna)
        object.__setattr__(self, "usco_hit_count", usco_hit_count)
        object.__setattr__(self, "usco_genes", usco_genes)
        object.__setattr__(self, "tiara_class", tiara_class)
        object.__setattr__(self, "telomere_offset_bp", telomere_offset_bp)
        object.__setattr__(self, "telomere_tract_bp", telomere_tract_bp)

        is_5p = end == "5p" or (port.strand == Strand.MINUS)
        if telomere_canonical_count is not None and telomere_opposite_count is not None:
            canon = telomere_canonical_count
            opp = telomere_opposite_count
            fwd = opp if is_5p else canon
            rev = canon if is_5p else opp
            if telomere_fwd_count is not None:
                fwd = telomere_fwd_count
            if telomere_rev_count is not None:
                rev = telomere_rev_count
        else:
            fwd = telomere_fwd_count if telomere_fwd_count is not None else 0
            rev = telomere_rev_count if telomere_rev_count is not None else 0
            canon = rev if is_5p else fwd
            opp = fwd if is_5p else rev

        object.__setattr__(self, "telomere_canonical_count", canon)
        object.__setattr__(self, "telomere_opposite_count", opp)
        object.__setattr__(self, "telomere_fwd_count", fwd)
        object.__setattr__(self, "telomere_rev_count", rev)


def build_end_evidence_ledger(
    graph: BidirectedGraph,
    annotations: dict[str, NodeAnnotations],
    terminal_window_bp: int | None = None,
    policy: Policy | None = None,
) -> list[EndEvidence]:
    """Build end-evidence ledger evaluating both ends (5p and 3p) of all segments in graph.

    Args:
        graph: Bidirected assembly graph.
        annotations: Node annotations mapping sequence_id to NodeAnnotations.
        terminal_window_bp: Terminal window size in bp. If None, uses default.
        policy: Policy instance governing thresholds. If None, uses DEFAULT_POLICY.

    Returns:
        Deterministically sorted list of EndEvidence records by (segment_name, end).
    """
    active_policy = policy if policy is not None else DEFAULT_POLICY
    min_tel_repeats = active_policy.telomere.min_terminal_repeat_count.value
    window_bp = terminal_window_bp if terminal_window_bp is not None else DEFAULT_TERMINAL_WINDOW_BP
    ledger: list[EndEvidence] = []

    for seg_name, seg in graph.segments.items():
        ann = annotations.get(seg_name)
        seg_len = seg.length
        seg_depth = seg.depth

        is_short = seg_len < 2 * window_bp
        midpoint = seg_len // 2

        boundary_5p = midpoint if is_short else window_bp
        boundary_3p = midpoint if is_short else max(0, seg_len - window_bp)

        if ann is not None:
            usco_count = len(ann.compleasm_hits)
            usco_genes = tuple(sorted({h.gene_id for h in ann.compleasm_hits}))
            tiara_cls = ann.tiara_classification if ann.tiara_classification is not None else ""
        else:
            usco_count = 0
            usco_genes = ()
            tiara_cls = ""

        # 5' port: SegmentPort(segment.name, Strand.MINUS)
        p_5p = SegmentPort(seg_name, Strand.MINUS)
        out_5p = graph.get_outgoing(p_5p)
        deg_5p = len(out_5p)
        tip_5p = deg_5p == 0
        neigh_5p = tuple(dict.fromkeys(e.to_port.to_symbol() for e in out_5p))

        # 3' port: SegmentPort(segment.name, Strand.PLUS)
        p_3p = SegmentPort(seg_name, Strand.PLUS)
        out_3p = graph.get_outgoing(p_3p)
        deg_3p = len(out_3p)
        tip_3p = deg_3p == 0
        neigh_3p = tuple(dict.fromkeys(e.to_port.to_symbol() for e in out_3p))

        seq = graph.get_sequence(seg_name)
        if seq:
            max_offset = active_policy.telomere.max_terminal_offset_bp.value
            max_gap = active_policy.telomere.max_internal_gap_bp.value
            scan_bp = boundary_5p

            call_5p_canon = scan_terminal_tract(
                sequence=seq,
                motif=active_policy.telomere.motif_reverse.value,
                terminal_scan_bp=scan_bp,
                max_terminal_offset_bp=max_offset,
                max_internal_gap_bp=max_gap,
                min_repeat_count=min_tel_repeats,
                is_suffix=False,
            )
            call_5p_opp = scan_terminal_tract(
                sequence=seq,
                motif=active_policy.telomere.motif_forward.value,
                terminal_scan_bp=scan_bp,
                max_terminal_offset_bp=max_offset,
                max_internal_gap_bp=max_gap,
                min_repeat_count=min_tel_repeats,
                is_suffix=False,
            )
            tel_canon_5p = call_5p_canon.repeat_count
            tel_opp_5p = call_5p_opp.repeat_count
            has_tel_5p = call_5p_canon.is_terminal_cap
            tel_fwd_5p = tel_opp_5p
            tel_rev_5p = tel_canon_5p
            tel_offset_5p = call_5p_canon.offset_bp
            tel_tract_5p = call_5p_canon.tract_length_bp

            call_3p_canon = scan_terminal_tract(
                sequence=seq,
                motif=active_policy.telomere.motif_forward.value,
                terminal_scan_bp=scan_bp,
                max_terminal_offset_bp=max_offset,
                max_internal_gap_bp=max_gap,
                min_repeat_count=min_tel_repeats,
                is_suffix=True,
            )
            call_3p_opp = scan_terminal_tract(
                sequence=seq,
                motif=active_policy.telomere.motif_reverse.value,
                terminal_scan_bp=scan_bp,
                max_terminal_offset_bp=max_offset,
                max_internal_gap_bp=max_gap,
                min_repeat_count=min_tel_repeats,
                is_suffix=True,
            )
            tel_canon_3p = call_3p_canon.repeat_count
            tel_opp_3p = call_3p_opp.repeat_count
            has_tel_3p = call_3p_canon.is_terminal_cap
            tel_fwd_3p = tel_canon_3p
            tel_rev_3p = tel_opp_3p
            tel_offset_3p = call_3p_canon.offset_bp
            tel_tract_3p = call_3p_canon.tract_length_bp
        else:
            tel_canon_5p = 0
            tel_opp_5p = 0
            tel_fwd_5p = 0
            tel_rev_5p = 0
            has_tel_5p = False
            tel_offset_5p = None
            tel_tract_5p = None

            tel_canon_3p = 0
            tel_opp_3p = 0
            tel_fwd_3p = 0
            tel_rev_3p = 0
            has_tel_3p = False
            tel_offset_3p = None
            tel_tract_3p = None

        if ann is not None:
            rdna_5p = [h for h in ann.barrnap_hits if h.start <= boundary_5p]
            rdna_cnt_5p = len(rdna_5p)
            rdna_feat_5p = tuple(sorted({h.feature_name for h in rdna_5p}))
            has_rdna_5p = rdna_cnt_5p > 0

            rdna_3p = [h for h in ann.barrnap_hits if h.end >= boundary_3p]
            rdna_cnt_3p = len(rdna_3p)
            rdna_feat_3p = tuple(sorted({h.feature_name for h in rdna_3p}))
            has_rdna_3p = rdna_cnt_3p > 0
        else:
            rdna_cnt_5p = 0
            rdna_feat_5p = ()
            has_rdna_5p = False
            rdna_cnt_3p = 0
            rdna_feat_3p = ()
            has_rdna_3p = False

        ev_5p = EndEvidence(
            segment_name=seg_name,
            end="5p",
            port=p_5p,
            length=seg_len,
            depth=seg_depth,
            degree=deg_5p,
            is_tip=tip_5p,
            outgoing_neighbors=neigh_5p,
            telomere_canonical_count=tel_canon_5p,
            telomere_opposite_count=tel_opp_5p,
            has_terminal_telomere=has_tel_5p,
            rdna_hit_count=rdna_cnt_5p,
            rdna_features=rdna_feat_5p,
            has_terminal_rdna=has_rdna_5p,
            usco_hit_count=usco_count,
            usco_genes=usco_genes,
            tiara_class=tiara_cls,
            telomere_fwd_count=tel_fwd_5p,
            telomere_rev_count=tel_rev_5p,
            telomere_offset_bp=tel_offset_5p,
            telomere_tract_bp=tel_tract_5p,
        )
        ledger.append(ev_5p)

        ev_3p = EndEvidence(
            segment_name=seg_name,
            end="3p",
            port=p_3p,
            length=seg_len,
            depth=seg_depth,
            degree=deg_3p,
            is_tip=tip_3p,
            outgoing_neighbors=neigh_3p,
            telomere_canonical_count=tel_canon_3p,
            telomere_opposite_count=tel_opp_3p,
            has_terminal_telomere=has_tel_3p,
            rdna_hit_count=rdna_cnt_3p,
            rdna_features=rdna_feat_3p,
            has_terminal_rdna=has_rdna_3p,
            usco_hit_count=usco_count,
            usco_genes=usco_genes,
            tiara_class=tiara_cls,
            telomere_fwd_count=tel_fwd_3p,
            telomere_rev_count=tel_rev_3p,
            telomere_offset_bp=tel_offset_3p,
            telomere_tract_bp=tel_tract_3p,
        )
        ledger.append(ev_3p)

    ledger.sort(key=lambda e: (e.segment_name, e.end))
    return ledger


def emit_tsv_ledger(ledger: list[EndEvidence], out_file: TextIO | Path | str) -> None:
    """Emit standardized TSV ledger with header.

    Args:
        ledger: List of EndEvidence records to serialize.
        out_file: File path or writable text stream.
    """
    header = [
        "segment_name",
        "end",
        "port",
        "length",
        "depth",
        "degree",
        "is_tip",
        "outgoing_neighbors",
        "telomere_canonical_count",
        "telomere_opposite_count",
        "telomere_fwd_count",
        "telomere_rev_count",
        "has_terminal_telomere",
        "rdna_hit_count",
        "rdna_features",
        "has_terminal_rdna",
        "usco_hit_count",
        "usco_genes",
        "tiara_class",
    ]

    def _write_content(writer: Any) -> None:
        writer.write("\t".join(header) + "\n")
        for e in ledger:
            row = [
                e.segment_name,
                e.end,
                e.port.to_symbol(),
                str(e.length),
                str(e.depth),
                str(e.degree),
                str(e.is_tip),
                ",".join(e.outgoing_neighbors),
                str(e.telomere_canonical_count),
                str(e.telomere_opposite_count),
                str(e.telomere_fwd_count),
                str(e.telomere_rev_count),
                str(e.has_terminal_telomere),
                str(e.rdna_hit_count),
                ",".join(e.rdna_features),
                str(e.has_terminal_rdna),
                str(e.usco_hit_count),
                ",".join(e.usco_genes),
                e.tiara_class,
            ]
            writer.write("\t".join(row) + "\n")

    if isinstance(out_file, str | Path):
        with open(out_file, "w", encoding="utf-8") as f:
            _write_content(f)
    else:
        _write_content(out_file)
