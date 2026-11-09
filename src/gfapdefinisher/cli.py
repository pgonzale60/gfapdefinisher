"""Command-line interface for gfapdefinisher."""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from gfapdefinisher.annotations import aggregate_node_annotations
from gfapdefinisher.cigar import parse_dovetail_overlap_len
from gfapdefinisher.classifier import _lookup_evidence, compute_core_coverage
from gfapdefinisher.emitter import (
    DEFAULT_AGP_POLICY,
    DEFAULT_MIN_SUPER_LENGTH,
    AgpRecord,
    PathComponent,
    SequenceProvenance,
    assemble_scaffold_sequence,
    compute_break_site_summary,
    emit_agp,
    emit_excluded_fasta,
    emit_fasta,
    emit_grs_bed,
    emit_summary_tsv,
    emit_telomere_bed,
    format_summary_table,
    nominate_super_scaffolds,
)
from gfapdefinisher.graph import BidirectedGraph, load_gfa
from gfapdefinisher.ledger import EndEvidence, build_end_evidence_ledger, emit_tsv_ledger
from gfapdefinisher.models import SegmentPort, Strand
from gfapdefinisher.policy import DEFAULT_POLICY
from gfapdefinisher.rdna import RdnaCassetteResolver, RdnaResolution
from gfapdefinisher.resolutions import (
    EndState,
    EndStateRecord,
    JunctionRecord,
    ResolutionRecord,
    ResolutionType,
    build_feature_id,
    collect_features,
    emit_end_states_tsv,
    emit_junctions_tsv,
    emit_resolutions_tsv,
)
from gfapdefinisher.resolver import (
    OntBridgingEngine,
    SatelliteGapResolution,
    _parse_port_spec,
    reverse_complement,
)


def create_parser() -> argparse.ArgumentParser:
    """Create and configure top-level argument parser."""
    parser = argparse.ArgumentParser(
        prog="gfapdefinisher",
        description="Biologically informed Hifiasm GFA harness for automated germline T2T assemblies",
    )
    subparsers = parser.add_subparsers(dest="command", help="Subcommand to execute")

    report_parser = subparsers.add_parser("report", help="Generate unitig end-evidence ledger TSV report")
    report_parser.add_argument("--gfa", required=True, type=Path, help="Path to input GFA1 file (.gfa or .gfa.gz)")
    report_parser.add_argument("--tiara", type=Path, default=None, help="Path to Tiara classification tabular output")
    report_parser.add_argument("--barrnap", type=Path, default=None, help="Path to Barrnap rDNA annotations GFF3")
    report_parser.add_argument("--out", type=Path, default=None, help="Output TSV path (default: stdout)")
    report_parser.add_argument(
        "--terminal-window",
        type=int,
        default=None,
        dest="terminal_window",
        help="Terminal window in bp to evaluate end-evidence",
    )

    finish_parser = subparsers.add_parser(
        "finish", help="End-to-end pipeline CLI for T2T finishing and ledger emission"
    )
    finish_parser.add_argument("--gfa", required=True, type=Path, help="Path to input GFA1 file (.gfa or .gfa.gz)")
    finish_parser.add_argument("--tiara", type=Path, default=None, help="Path to Tiara classification tabular output")
    finish_parser.add_argument("--barrnap", type=Path, default=None, help="Path to Barrnap rDNA annotations GFF3")
    finish_parser.add_argument(
        "--ont-paf",
        type=Path,
        default=None,
        dest="ont_paf",
        help="Optional path to spanning long-read PAF alignment",
    )
    finish_parser.add_argument(
        "--ont-reads",
        type=Path,
        default=None,
        dest="ont_reads",
        help="Optional path to long reads FASTA/FASTQ for representative sequence slicing",
    )
    finish_parser.add_argument(
        "--out-dir",
        required=True,
        type=Path,
        dest="out_dir",
        help="Directory where all ledger TSVs and artifacts will be written",
    )
    finish_parser.add_argument(
        "--terminal-window",
        type=int,
        default=None,
        dest="terminal_window",
        help="Window size in bp to evaluate end evidence",
    )
    finish_parser.add_argument(
        "--emit-fasta",
        nargs="?",
        const=True,
        default=None,
        dest="emit_fasta",
        help="Optional path to emit assembled scaffold FASTA (default: <out_dir>/scaffolds.fasta)",
    )
    finish_parser.add_argument(
        "--emit-agp",
        nargs="?",
        const=True,
        default=None,
        dest="emit_agp",
        help="Optional path to emit AGP v2.1 records (default: <out_dir>/scaffolds.agp)",
    )
    finish_parser.add_argument(
        "--excluded-fasta",
        nargs="?",
        const=True,
        default=None,
        dest="excluded_fasta",
        help="Optional path to emit excluded unitigs FASTA (default: <out_dir>/excluded.fasta)",
    )
    finish_parser.add_argument(
        "--resolutions",
        type=Path,
        default=None,
        dest="resolutions",
        help="Optional custom output path for resolutions TSV (default: <out_dir>/resolutions.tsv)",
    )
    finish_parser.add_argument(
        "--prefix",
        type=str,
        default="",
        help="Optional prefix for scaffold identifiers",
    )
    finish_parser.add_argument(
        "--nominate-super",
        action="store_true",
        default=False,
        dest="nominate_super",
        help="Nominate large sequences with break site evidence as SUPER chromosomes",
    )
    finish_parser.add_argument(
        "--min-super-length",
        "--min-super-len",
        type=int,
        default=DEFAULT_MIN_SUPER_LENGTH,
        dest="min_super_length",
        help="Minimum length in bp for a sequence with break site evidence to be nominated as a SUPER chromosome (default: 100000)",
    )
    finish_parser.add_argument(
        "--emit-grs-bed",
        nargs="?",
        const=True,
        default=None,
        dest="emit_grs_bed",
        help="Optional path to emit GRS regions BED (default: <out_dir>/GRS.bed)",
    )
    finish_parser.add_argument(
        "--emit-telomere-bed",
        nargs="?",
        const=True,
        default=None,
        dest="emit_telomere_bed",
        help="Optional path to emit telomere caps BED (default: <out_dir>/telomeres.bed)",
    )
    finish_parser.add_argument(
        "--summary",
        nargs="?",
        const=True,
        default=None,
        dest="summary",
        help="Optional path to emit break-site summary TSV (default: <out_dir>/summary.tsv)",
    )

    return parser


def _parse_port(p_str: str) -> SegmentPort:
    name, strand = _parse_port_spec(p_str)
    return SegmentPort(name, strand)


def _find_overlap_len(graph: BidirectedGraph, from_port: SegmentPort, to_port: SegmentPort) -> int:
    """Find effective overlap length from from_port to to_port in dovetail edges."""
    edges = graph.edges_from(from_port) if hasattr(graph, "edges_from") else graph.get_outgoing(from_port)
    for e in edges:
        if e.to_port == to_port:
            return parse_dovetail_overlap_len(e.cigar) if (e.cigar and e.cigar != "*") else e.overlap_len
    return 0


def _is_contam_unitig(ledger_map: Mapping[str, Any], name: str) -> bool:
    contam_classes = ("bacteria", "archaea", "prokarya", "contaminant")
    for key in (f"{name}+", f"{name}-", f"{name}:5p", f"{name}:3p", name):
        ev = ledger_map.get(key)
        if ev is not None:
            if getattr(ev, "is_contaminant", False):
                return True
            cls_str = str(getattr(ev, "tiara_class", "")).lower().strip()
            if cls_str in contam_classes:
                return True
    return False


def _get_usco_count(ledger_map: Mapping[str, Any], name: str) -> int:
    for key in (f"{name}+", f"{name}-", f"{name}:5p", f"{name}:3p", name):
        ev = ledger_map.get(key)
        if ev is not None and getattr(ev, "usco_hit_count", 0) > 0:
            return ev.usco_hit_count
    return 0


def _trace_germline_chain(
    graph: BidirectedGraph,
    ledger_map: Mapping[str, Any],
    germ_port: SegmentPort,
    max_hops: int,
    core_coverage: float = 0.0,
) -> list[SegmentPort]:
    """Trace outward germline continuation unitigs to an authentic terminal telomeric cap."""
    ev = _lookup_evidence(ledger_map, germ_port)
    deg = graph.degree(germ_port) if germ_port.name in graph.segments else (ev.degree if ev else 0)
    if ev is not None and ev.has_terminal_telomere and deg == 0:
        return [germ_port]

    if germ_port.name not in graph.segments:
        return [germ_port]

    queue: list[tuple[SegmentPort, list[SegmentPort]]] = [(germ_port, [germ_port])]
    visited: set[SegmentPort] = {germ_port}

    min_core_ratio = DEFAULT_POLICY.stoichiometry.min_stem_core_ratio.value

    while queue:
        curr, path = queue.pop(0)
        if len(path) - 1 >= max_hops:
            continue

        edges = graph.edges_from(curr) if hasattr(graph, "edges_from") else graph.get_outgoing(curr)
        for e in edges:
            nxt = e.to_port
            if nxt in visited:
                continue
            nxt_seg = graph.segments.get(nxt.name)
            is_core_depth = (
                (nxt_seg.depth or 0.0) >= min_core_ratio * core_coverage
                if nxt_seg is not None and core_coverage > 0
                else False
            )
            if is_core_depth or _get_usco_count(ledger_map, nxt.name) != 0 or _is_contam_unitig(ledger_map, nxt.name):
                continue

            nxt_ev = _lookup_evidence(ledger_map, nxt)
            nxt_deg = graph.degree(nxt) if nxt.name in graph.segments else (nxt_ev.degree if nxt_ev else 0)
            if nxt_ev is not None and nxt_ev.has_terminal_telomere and nxt_deg == 0:
                return [*path, nxt]

            visited.add(nxt)
            queue.append((nxt, [*path, nxt]))

    return [germ_port]


def build_scaffolds(
    graph: BidirectedGraph,
    ledger_map: Mapping[str, Any],
    end_states: Sequence[EndStateRecord],
    resolutions: Sequence[ResolutionRecord] | None = None,
    gap_resolutions: Sequence[SatelliteGapResolution] = (),
    dropout_resolutions: Sequence[tuple[ResolutionRecord, EndStateRecord, str | None]] = (),
    rdna_resolutions: Sequence[RdnaResolution] = (),
    prefix: str = "",
    junctions: Sequence[JunctionRecord] = (),
    *,
    quirks: Sequence[ResolutionRecord] | None = None,
    resolution_records: list[ResolutionRecord] | None = None,
    nominate_super: bool = False,
    min_super_length: int = DEFAULT_MIN_SUPER_LENGTH,
    core_coverage: float = 0.0,
) -> tuple[dict[str, str], list[AgpRecord]]:
    """Assemble scaffold sequences and AGP v2.1 records with empirical provenance."""
    if core_coverage <= 0.0 and ledger_map:
        ev_items = [v for v in ledger_map.values() if isinstance(v, EndEvidence)]
        if ev_items:
            core_coverage = compute_core_coverage(graph, ev_items)

    active_resolutions: list[ResolutionRecord] | None = None
    if resolution_records is not None:
        active_resolutions = resolution_records
    elif isinstance(resolutions, list):
        active_resolutions = resolutions
    elif isinstance(quirks, list):
        active_resolutions = quirks

    rdna_5p: dict[str, RdnaResolution] = {}
    rdna_3p: dict[str, RdnaResolution] = {}
    for r in rdna_resolutions:
        if len(r.path) > 1:
            if r.side == "5p":
                rdna_5p[r.path[0].name] = r
            else:
                rdna_3p[r.path[0].name] = r

    gap_by_left: dict[tuple[str, str], SatelliteGapResolution] = {}
    gap_by_left_seg: dict[str, SatelliteGapResolution] = {}
    for g in gap_resolutions:
        r_rec = g.resolution_record
        if r_rec and r_rec.left_anchor_port:
            l_name, l_strand = _parse_port_spec(r_rec.left_anchor_port)
            gap_by_left[(l_name, l_strand.value)] = g
            gap_by_left_seg[l_name] = g
        elif g.chrom:
            gap_by_left_seg[g.chrom] = g

    dropout_by_anchor: dict[str, tuple[ResolutionRecord, EndStateRecord, str]] = {}
    for rrec, erec, dseq in dropout_resolutions:
        if dseq:
            dropout_by_anchor[erec.chrom] = (rrec, erec, dseq)

    junctions_by_stem: dict[str, JunctionRecord] = {j.stem_port: j for j in junctions}
    junctions_by_germ: dict[str, JunctionRecord] = {
        j.germline_continuation_port: j for j in junctions if j.germline_continuation_port
    }
    max_hops = DEFAULT_POLICY.pde.max_somatic_cap_hops.value

    # Identify intercalary bridge linkages between stems
    bridge_forward: dict[SegmentPort, tuple[SegmentPort, SegmentPort]] = {}
    has_incoming_bridge: set[SegmentPort] = set()

    stem_names = {j.stem_port.rstrip("+-") for j in junctions}

    min_core_ratio = DEFAULT_POLICY.stoichiometry.min_stem_core_ratio.value

    def _has_telomere(port: SegmentPort) -> bool:
        ev = ledger_map.get(port.to_symbol())
        return bool(ev and getattr(ev, "has_terminal_telomere", False))

    def _is_core(name: str) -> bool:
        seg = graph.segments.get(name)
        is_core_depth = (
            (seg.depth or 0.0) >= min_core_ratio * core_coverage if seg is not None and core_coverage > 0 else False
        )
        return (name in stem_names or is_core_depth or _get_usco_count(ledger_map, name) > 0) and not _is_contam_unitig(
            ledger_map, name
        )

    seen_pde_stems = {
        (rec.left_anchor_port.rstrip("+-"), rec.side)
        for rec in (active_resolutions or [])
        if rec.resolution_type == ResolutionType.PDE_TERMINAL_CAP and rec.left_anchor_port
    } | {
        (rec.chrom, rec.side)
        for rec in (active_resolutions or [])
        if rec.resolution_type == ResolutionType.PDE_TERMINAL_CAP
    }

    seen_bridges: set[str] = set()
    for j in junctions:
        if j.germline_continuation_port:
            stem_p = _parse_port(j.stem_port)
            gp = _parse_port(j.germline_continuation_port)
            gp_inv_sym = gp.invert().to_symbol()
            if gp_inv_sym in junctions_by_germ:
                j2 = junctions_by_germ[gp_inv_sym]
                next_stem_p = _parse_port(j2.stem_port).invert()
                bridge_forward[stem_p] = (gp, next_stem_p)
                has_incoming_bridge.add(next_stem_p)
                for j_elem in (j, j2):
                    seen_pde_stems.add((j_elem.chrom, j_elem.side))
                    seen_pde_stems.add((j_elem.stem_port.rstrip("+-"), j_elem.side))

                if active_resolutions is not None and gp.name not in seen_bridges:
                    seen_bridges.add(gp.name)
                    bridge_count = len(
                        [
                            rec
                            for rec in active_resolutions
                            if rec.resolution_type == ResolutionType.INTERCALARY_GRS_BRIDGE
                        ]
                    )
                    bridge_fid = build_feature_id(
                        ResolutionType.INTERCALARY_GRS_BRIDGE,
                        stem_p.name,
                        "mid",
                        index=bridge_count + 1,
                    )
                    som_j = j.somatic_terminus_port.rstrip("+-") if j.somatic_terminus_port else None
                    som_j2 = j2.somatic_terminus_port.rstrip("+-") if j2.somatic_terminus_port else None
                    excl_stubs = tuple(s for s in (som_j, som_j2) if s)
                    bridge_rec = ResolutionRecord(
                        feature_id=bridge_fid,
                        chrom=stem_p.name,
                        side="mid",
                        resolution_type=ResolutionType.INTERCALARY_GRS_BRIDGE,
                        left_anchor_port=j.stem_port,
                        right_anchor_port=j2.stem_port,
                        unitig_ids=(gp.name,),
                        excluded_unitig_ids=excl_stubs,
                        core_coverage_bp_median=(graph.segments[gp.name].depth if gp.name in graph.segments else None),
                        notes=(
                            f"Intercalary GRS bridge {gp.to_symbol()} connecting core stems "
                            f"{j.stem_port} and {j2.stem_port}"
                        ),
                    )
                    active_resolutions.append(bridge_rec)

    if active_resolutions is not None:
        existing_pde_caps: dict[tuple[str, str], ResolutionRecord] = {}
        for rec in active_resolutions:
            if rec.resolution_type == ResolutionType.PDE_TERMINAL_CAP:
                if rec.left_anchor_port:
                    existing_pde_caps[(rec.left_anchor_port.rstrip("+-"), rec.side)] = rec
                existing_pde_caps[(rec.chrom, rec.side)] = rec

        for j in junctions:
            if not j.germline_continuation_port:
                continue
            gp = _parse_port(j.germline_continuation_port)
            if gp.invert().to_symbol() in junctions_by_germ:
                continue
            cap_path = _trace_germline_chain(graph, ledger_map, gp, max_hops, core_coverage=core_coverage)
            has_telo = bool(cap_path and _has_telomere(cap_path[-1]))
            target_ports = cap_path if has_telo else [gp]
            som_excl = (j.somatic_terminus_port.rstrip("+-"),) if j.somatic_terminus_port else ()

            j_key = (j.chrom, j.side)
            stem_key = (j.stem_port.rstrip("+-"), j.side)
            ex_rec = existing_pde_caps.get(j_key) or existing_pde_caps.get(stem_key)

            if ex_rec is not None:
                # Update existing record with the full multi-hop chain
                ex_rec.unitig_ids = tuple(p.name for p in target_ports)
                ex_rec.right_anchor_port = target_ports[-1].to_symbol()
                if som_excl and not ex_rec.excluded_unitig_ids:
                    ex_rec.excluded_unitig_ids = som_excl
            elif j_key not in seen_pde_stems and stem_key not in seen_pde_stems:
                seen_pde_stems.add(j_key)
                seen_pde_stems.add(stem_key)
                j_fid = build_feature_id(ResolutionType.PDE_TERMINAL_CAP, j.chrom, j.side, index=1)
                j_rec = ResolutionRecord(
                    feature_id=j_fid,
                    chrom=j.chrom,
                    side=j.side,
                    resolution_type=ResolutionType.PDE_TERMINAL_CAP,
                    left_anchor_port=j.stem_port,
                    right_anchor_port=target_ports[-1].to_symbol(),
                    unitig_ids=tuple(p.name for p in target_ports),
                    excluded_unitig_ids=som_excl,
                    core_coverage_bp_median=(graph.segments[j.chrom].depth if j.chrom in graph.segments else None),
                    notes=(
                        f"Resolved terminal PDE junction via germline continuation chain "
                        f"{' -> '.join(p.to_symbol() for p in target_ports)} to telomeric cap "
                        f"{target_ports[-1].to_symbol()}; dismissed somatic sister stub "
                        f"{j.somatic_terminus_port or ''}"
                        if has_telo and len(target_ports) > 1
                        else f"Resolved terminal PDE junction via germline continuation {gp.to_symbol()}; "
                        f"dismissed somatic sister stub {j.somatic_terminus_port or ''}"
                    ),
                )
                active_resolutions.append(j_rec)
                existing_pde_caps[j_key] = j_rec
                existing_pde_caps[stem_key] = j_rec

    consumed: set[str] = set()
    if active_resolutions is not None:
        for rec in active_resolutions:
            consumed.update(rec.excluded_unitig_ids)
    elif quirks is not None:
        for rec in quirks:
            consumed.update(rec.excluded_unitig_ids)
    for r in rdna_resolutions:
        consumed.update(r.excluded_unitigs)

    continuation_targets: set[str] = set()
    for gres in gap_resolutions:
        for ip in gres.intermediate_ports:
            continuation_targets.add(ip.name)
        if gres.right_flank_port:
            continuation_targets.add(gres.right_flank_port.name)
        elif gres.resolution_record and gres.resolution_record.right_anchor_port:
            r_name = gres.resolution_record.right_anchor_port.rstrip("+-")
            continuation_targets.add(r_name)

    for j in junctions:
        if j.germline_continuation_port:
            gp = _parse_port(j.germline_continuation_port)
            g_chain = _trace_germline_chain(graph, ledger_map, gp, max_hops, core_coverage=core_coverage)
            for p in g_chain:
                continuation_targets.add(p.name)

    for _stem_p, (bridge_p, next_stem) in bridge_forward.items():
        continuation_targets.add(bridge_p.name)
        continuation_targets.add(next_stem.name)

    for _rrec, erec, _dseq in dropout_resolutions:
        if erec.germline_port:
            g_name = erec.germline_port.rstrip("+-")
            continuation_targets.add(g_name)

    for r in rdna_resolutions:
        if len(r.path) > 1:
            for p in r.path[1:]:
                continuation_targets.add(p.name)
                consumed.add(p.name)
        if r.distal_cap_port is not None:
            continuation_targets.add(r.distal_cap_port.name)
            consumed.add(r.distal_cap_port.name)
        if r.is_intercalary and r.exit_stem_port is not None:
            continuation_targets.add(r.exit_stem_port.name)

    continuation_targets: set[str] = {s for s in continuation_targets if not _is_core(s)}

    seen_bubbles: set[tuple[str, str]] = set()

    def _step_forward(curr: SegmentPort, visited_stems: set[str]) -> list[SegmentPort] | None:
        if curr in bridge_forward:
            bridge_p, next_stem = bridge_forward[curr]
            if next_stem.name not in visited_stems and next_stem.name not in consumed and bridge_p.name not in consumed:
                return [bridge_p, next_stem]
        outgoing = [
            e.to_port
            for e in graph.get_outgoing(curr)
            if not _is_contam_unitig(ledger_map, e.to_port.name)
            and e.to_port.name not in visited_stems
            and e.to_port.name not in consumed
        ]
        cand_branches = [p for p in outgoing if not _is_core(p.name) and not _has_telomere(p)]
        if 0 < len(cand_branches) <= 2 and len(cand_branches) == len(outgoing):
            if len(cand_branches) == 1:
                b = cand_branches[0]
                b_outs = [
                    e.to_port
                    for e in graph.get_outgoing(b)
                    if not _is_contam_unitig(ledger_map, e.to_port.name)
                    and e.to_port.name not in visited_stems
                    and e.to_port.name not in consumed
                ]
                core_outs = [p for p in b_outs if _is_core(p.name)]
                if len(core_outs) == 1:
                    return [b, core_outs[0]]
            elif len(cand_branches) == 2:
                reconv_ports = set()
                for b in cand_branches:
                    b_outs = [
                        e.to_port
                        for e in graph.get_outgoing(b)
                        if not _is_contam_unitig(ledger_map, e.to_port.name)
                        and e.to_port.name not in visited_stems
                        and e.to_port.name not in consumed
                    ]
                    if len(b_outs) == 1:
                        reconv_ports.add(b_outs[0])
                if len(reconv_ports) == 1:
                    reconv_port = next(iter(reconv_ports))
                    if (
                        not _has_telomere(reconv_port)
                        and not _is_contam_unitig(ledger_map, reconv_port.name)
                        and reconv_port.name not in visited_stems
                        and reconv_port.name not in consumed
                    ):
                        reconv_outs = [
                            e.to_port
                            for e in graph.get_outgoing(reconv_port)
                            if not _is_contam_unitig(ledger_map, e.to_port.name)
                            and e.to_port.name not in visited_stems
                            and e.to_port.name not in consumed
                        ]
                        core_outs = [p for p in reconv_outs if _is_core(p.name)]
                        if len(core_outs) == 1:
                            best_branch = max(
                                cand_branches,
                                key=lambda b: (
                                    graph.segments[b.name].depth or 0.0,
                                    graph.segments[b.name].length,
                                    b.name,
                                ),
                            )
                            for b in cand_branches:
                                consumed.add(b.name)
                            consumed.add(reconv_port.name)
                            bubble_key = (
                                min(curr.to_symbol(), core_outs[0].to_symbol()),
                                max(curr.to_symbol(), core_outs[0].to_symbol()),
                            )
                            if active_resolutions is not None and bubble_key not in seen_bubbles:
                                seen_bubbles.add(bubble_key)
                                excl_branches = tuple(b.name for b in cand_branches if b != best_branch)
                                for eb in excl_branches:
                                    consumed.add(eb)
                                bubble_count = len(
                                    [
                                        rec
                                        for rec in active_resolutions
                                        if rec.resolution_type == ResolutionType.SOMATIC_ALLELIC_BUBBLE
                                    ]
                                )
                                bubble_fid = build_feature_id(
                                    ResolutionType.SOMATIC_ALLELIC_BUBBLE,
                                    curr.name,
                                    "mid",
                                    index=bubble_count + 1,
                                )
                                bubble_rec = ResolutionRecord(
                                    feature_id=bubble_fid,
                                    chrom=curr.name,
                                    side="mid",
                                    resolution_type=ResolutionType.SOMATIC_ALLELIC_BUBBLE,
                                    left_anchor_port=curr.to_symbol(),
                                    right_anchor_port=core_outs[0].to_symbol(),
                                    unitig_ids=(best_branch.name, reconv_port.name),
                                    excluded_unitig_ids=excl_branches,
                                    core_coverage_bp_median=(
                                        graph.segments[best_branch.name].depth
                                        if best_branch.name in graph.segments
                                        else None
                                    ),
                                    notes=(
                                        f"Traversed 2-hop somatic allelic bubble from {curr.to_symbol()} to {core_outs[0].to_symbol()}; "
                                        f"selected primary higher-depth allele {best_branch.to_symbol()} "
                                        f"({graph.segments[best_branch.name].depth if best_branch.name in graph.segments else 0.0:.1f}x) "
                                        f"over {','.join(excl_branches)}"
                                    ),
                                )
                                active_resolutions.append(bubble_rec)
                            return [best_branch, reconv_port, core_outs[0]]
        return None

    chain_by_seg: dict[str, list[SegmentPort]] = {}
    bridge_seg_names: set[str] = set()
    seen_chain_segs: set[str] = set()

    all_stem_ports = sorted(
        [
            p
            for p in bridge_forward
            if p not in has_incoming_bridge
            and (p.name not in continuation_targets or _is_core(p.name))
            and p.name not in consumed
        ]
        + [
            SegmentPort(s, Strand.FORWARD)
            for s in sorted(graph.segments.keys(), key=lambda s: (-graph.segments[s].length, s))
            if _is_core(s)
        ],
        key=lambda p: (-graph.segments[p.name].length if p.name in graph.segments else 0, p.name),
    )

    for start_p in all_stem_ports:
        if start_p.name in seen_chain_segs:
            continue
        chain: list[SegmentPort] = [start_p]
        visited_stems: set[str] = {start_p.name}
        curr = start_p
        while True:
            steps = _step_forward(curr, visited_stems)
            if not steps:
                break
            for p in steps[:-1]:
                bridge_seg_names.add(p.name)
            next_stem = steps[-1]
            visited_stems.add(next_stem.name)
            chain.extend(steps)
            curr = next_stem

        # Also walk backward from start_p
        curr_bwd = start_p.invert()
        bwd_chain: list[SegmentPort] = []
        while True:
            steps = _step_forward(curr_bwd, visited_stems)
            if not steps:
                break
            for p in steps[:-1]:
                bridge_seg_names.add(p.name)
            prev_stem_inv = steps[-1]
            visited_stems.add(prev_stem_inv.name)
            bwd_chain = [p.invert() for p in reversed(steps)] + bwd_chain
            curr_bwd = prev_stem_inv

        chain = bwd_chain + chain

        if len(chain) > 1:
            if chain[0].strand != Strand.FORWARD and any(p.strand == Strand.FORWARD for p in chain):
                chain = [p.invert() for p in reversed(chain)]
            seg_set = {p.name for p in chain}
            seen_chain_segs.update(seg_set)
            for p in chain:
                chain_by_seg[p.name] = chain

    scaffolds: dict[str, str] = {}
    all_agp: list[AgpRecord] = []
    placed_unitigs_all: set[str] = set()

    sorted_segs = sorted(graph.segments.keys(), key=lambda s: (-graph.segments[s].length, s))
    pass1_segs = [s for s in sorted_segs if s not in consumed and (_is_core(s) or s not in continuation_targets)]
    pass2_segs = [s for s in sorted_segs if s not in consumed and (not _is_core(s) and s in continuation_targets)]

    for seg_name in pass1_segs + pass2_segs:
        if seg_name in consumed or seg_name in placed_unitigs_all:
            continue
        seg = graph.segments[seg_name]
        if not seg.sequence:
            continue

        comps: list[PathComponent] = []

        def _add_comp(c: PathComponent) -> bool:
            if not c.is_gap:
                if c.segment_id in placed_unitigs_all:
                    return False
                placed_unitigs_all.add(c.segment_id)
                consumed.add(c.segment_id)
            comps.append(c)  # noqa: B023
            return True

        if seg_name in chain_by_seg:
            core_path = chain_by_seg[seg_name]
        else:
            core_path = [SegmentPort(seg_name, Strand.FORWARD)]

        first_port = core_path[0]
        last_port = core_path[-1]

        # 5p rDNA resolution
        if first_port.name in rdna_5p:
            r = rdna_5p[first_port.name]
            twin_path = [p.invert() for p in reversed(r.path)]
            if r.capped_via_heuristic and r.distal_cap_port is not None:
                cap_port = r.distal_cap_port
                cap_strand = cap_port.invert().strand.to_symbol()
                cap_seg = graph.segments.get(cap_port.name)
                cap_seq = (
                    cap_seg.sequence
                    if (cap_strand == "+" and cap_seg is not None)
                    else reverse_complement(cap_seg.sequence if cap_seg is not None else "")
                )
                if _add_comp(
                    PathComponent(
                        segment_id=cap_port.name,
                        strand=cap_strand,
                        sequence=cap_seq,
                        overlap_to_next_bp=0,
                        sequence_provenance=SequenceProvenance.HIFI_GERMLINE,
                    )
                ):
                    _add_comp(
                        PathComponent(
                            segment_id="rdna_tangle_gap",
                            strand="+",
                            is_gap=True,
                            gap_length=DEFAULT_AGP_POLICY.unknown_gap_bp.value,
                            gap_type="repeat",
                            linkage="yes",
                            linkage_evidence="rna_tangle",
                            overlap_to_next_bp=0,
                            gap_component_type="U",
                        )
                    )

            for i, port in enumerate(twin_path[:-1]):
                p_seg = graph.segments[port.name]
                p_seq = p_seg.sequence if port.strand == Strand.FORWARD else reverse_complement(p_seg.sequence or "")
                nxt_port = twin_path[i + 1] if i < len(twin_path) - 2 else first_port
                has_edge = any(e.to_port == nxt_port for e in graph.get_outgoing(port))
                overlap = _find_overlap_len(graph, port, nxt_port) if has_edge else 0
                prov = SequenceProvenance.RDNA_CHAINED
                if i == 0 and not (r.capped_via_heuristic and r.distal_cap_port is not None):
                    prov = SequenceProvenance.HIFI_GERMLINE
                if _add_comp(
                    PathComponent(
                        segment_id=port.name,
                        strand=port.strand.to_symbol(),
                        sequence=p_seq,
                        overlap_to_next_bp=overlap,
                        sequence_provenance=prov,
                    )
                ):
                    if not has_edge:
                        _add_comp(
                            PathComponent(
                                segment_id="rdna_tangle_gap",
                                strand="+",
                                is_gap=True,
                                gap_length=DEFAULT_AGP_POLICY.unknown_gap_bp.value,
                                gap_type="repeat",
                                linkage="yes",
                                linkage_evidence="rna_tangle",
                                overlap_to_next_bp=0,
                                gap_component_type="U",
                            )
                        )
            consumed.update(r.excluded_unitigs)
        else:
            # 5p germline continuation extension
            first_inv = first_port.invert()
            g_chain_5p: list[SegmentPort] = []
            if first_inv.to_symbol() in junctions_by_stem:
                j_5p = junctions_by_stem[first_inv.to_symbol()]
                if j_5p.germline_continuation_port:
                    gp_5p = _parse_port(j_5p.germline_continuation_port)
                    g_chain_5p = _trace_germline_chain(graph, ledger_map, gp_5p, max_hops, core_coverage=core_coverage)

            # 5p ONT gap resolution (multi-hop and unpatched satellite gaps)
            curr_out = g_chain_5p[-1] if g_chain_5p else first_port.invert()
            gaps_5p: list[SatelliteGapResolution] = []
            visited_5p_gaps: set[tuple[str, str]] = set()
            while True:
                key = (curr_out.name, curr_out.strand.value)
                if key in visited_5p_gaps:
                    break
                gres = gap_by_left.get(key)
                if gres is None:
                    cand = gap_by_left_seg.get(curr_out.name)
                    if cand is not None and cand.resolution_record and cand.resolution_record.side == "5p":
                        gres = cand
                if gres is None:
                    break
                visited_5p_gaps.add(key)
                gaps_5p.append(gres)
                if gres.right_flank_port:
                    curr_out = gres.right_flank_port
                elif gres.resolution_record and gres.resolution_record.right_anchor_port:
                    r_name, r_strand = _parse_port_spec(gres.resolution_record.right_anchor_port)
                    curr_out = SegmentPort(r_name, r_strand)
                else:
                    break

            for gres in reversed(gaps_5p):
                qrec = gres.resolution_record
                r_name = (
                    qrec.right_anchor_port.rstrip("+-")
                    if (qrec and qrec.right_anchor_port)
                    else (gres.right_flank_port.name if gres.right_flank_port else None)
                )
                if qrec is not None:
                    placed_u: list[str] = []
                    if r_name and r_name in graph.segments:
                        placed_u.append(r_name)
                    if gres.intermediate_ports:
                        for ip in reversed(gres.intermediate_ports):
                            if ip.name in graph.segments:
                                placed_u.append(ip.name)
                    elif gres.representative_sequence and gres.representative_read:
                        placed_u.append(gres.representative_read.qname)
                    if placed_u:
                        qrec.unitig_ids = tuple(placed_u)
                r_strand = (
                    Strand.FORWARD
                    if (qrec and qrec.right_anchor_port and "+" in qrec.right_anchor_port)
                    else (gres.right_flank_port.strand if gres.right_flank_port else Strand.FORWARD)
                )
                r_strand_inv = r_strand.invert()
                if r_name and r_name in graph.segments:
                    r_seg = graph.segments[r_name]
                    r_seq = (
                        r_seg.sequence if r_strand_inv == Strand.FORWARD else reverse_complement(r_seg.sequence or "")
                    )
                    _add_comp(
                        PathComponent(
                            segment_id=r_name,
                            strand=r_strand_inv.to_symbol(),
                            sequence=r_seq,
                            overlap_to_next_bp=0,
                            sequence_provenance=SequenceProvenance.HIFI_GERMLINE,
                        )
                    )

                if gres.intermediate_ports:
                    gap_ev = qrec.agp_linkage_evidence if (qrec and qrec.agp_linkage_evidence) else "satellite-kmer"
                    _add_comp(
                        PathComponent(
                            segment_id=f"gap_{gres.chrom}_{gres.side}_{gres.gap_index}_cap",
                            strand="+",
                            is_gap=True,
                            gap_length=DEFAULT_AGP_POLICY.unknown_gap_bp.value,
                            gap_type="repeat",
                            linkage="yes",
                            linkage_evidence=gap_ev,
                            overlap_to_next_bp=0,
                            gap_component_type="U",
                        )
                    )
                    for ip in reversed(gres.intermediate_ports):
                        if ip.name in graph.segments:
                            ip_seg = graph.segments[ip.name]
                            ip_seq = (
                                ip_seg.sequence
                                if ip.strand == Strand.FORWARD
                                else reverse_complement(ip_seg.sequence or "")
                            )
                            _add_comp(
                                PathComponent(
                                    segment_id=ip.name,
                                    strand=ip.strand.to_symbol(),
                                    sequence=ip_seq,
                                    overlap_to_next_bp=0,
                                    sequence_provenance=SequenceProvenance.HIFI_GERMLINE,
                                )
                            )
                            _add_comp(
                                PathComponent(
                                    segment_id=f"gap_{gres.chrom}_{gres.side}_{gres.gap_index}_anchor",
                                    strand="+",
                                    is_gap=True,
                                    gap_length=DEFAULT_AGP_POLICY.unknown_gap_bp.value,
                                    gap_type="repeat",
                                    linkage="yes",
                                    linkage_evidence=gap_ev,
                                    overlap_to_next_bp=0,
                                    gap_component_type="U",
                                )
                            )
                elif gres.representative_sequence and gres.representative_read:
                    _add_comp(
                        PathComponent(
                            segment_id=gres.representative_read.qname,
                            strand="-",
                            sequence=reverse_complement(gres.representative_sequence or ""),
                            overlap_to_next_bp=0,
                            sequence_provenance=SequenceProvenance.ONT_DERIVED,
                            unpolished=True,
                        )
                    )
                else:
                    gap_len = (
                        int(round(gres.median_gap_bp))
                        if gres.median_gap_bp
                        else DEFAULT_AGP_POLICY.unknown_gap_bp.value
                    )
                    gap_type = (
                        "repeat"
                        if (
                            qrec
                            and qrec.resolution_type
                            in (ResolutionType.SAT_GAP_VARIABLE, ResolutionType.SUBTELOMERIC_SATELLITE_CAP)
                        )
                        else "scaffold"
                    )
                    gap_ev = (
                        qrec.agp_linkage_evidence
                        if (qrec and qrec.agp_linkage_evidence)
                        else ("satellite-kmer" if gap_type == "repeat" else "hifi-spanning")
                    )
                    comp_type = "U" if gap_len == DEFAULT_AGP_POLICY.unknown_gap_bp.value else "N"
                    _add_comp(
                        PathComponent(
                            segment_id=f"gap_{gres.chrom}_{gres.side}_{gres.gap_index}",
                            strand="+",
                            is_gap=True,
                            gap_length=gap_len,
                            gap_type=gap_type,
                            linkage="yes",
                            linkage_evidence=gap_ev,
                            overlap_to_next_bp=0,
                            gap_component_type=comp_type,
                        )
                    )

            prefix_ports = [p.invert() for p in reversed(g_chain_5p)]
            for i, port in enumerate(prefix_ports):
                nxt_port = prefix_ports[i + 1] if i < len(prefix_ports) - 1 else first_port
                overlap = _find_overlap_len(graph, port, nxt_port)
                p_seg = graph.segments[port.name]
                p_seq = p_seg.sequence if port.strand == Strand.FORWARD else reverse_complement(p_seg.sequence or "")
                _add_comp(
                    PathComponent(
                        segment_id=port.name,
                        strand=port.strand.to_symbol(),
                        sequence=p_seq,
                        overlap_to_next_bp=overlap,
                        sequence_provenance=SequenceProvenance.HIFI_GERMLINE,
                    )
                )

        # Core path
        for i, port in enumerate(core_path):
            overlap = _find_overlap_len(graph, port, core_path[i + 1]) if i < len(core_path) - 1 else 0
            p_seg = graph.segments[port.name]
            p_seq = p_seg.sequence if port.strand == Strand.FORWARD else reverse_complement(p_seg.sequence or "")
            prov = SequenceProvenance.HIFI_GERMLINE if port.name in bridge_seg_names else SequenceProvenance.HIFI_CORE
            _add_comp(
                PathComponent(
                    segment_id=port.name,
                    strand=port.strand.to_symbol(),
                    sequence=p_seq,
                    overlap_to_next_bp=overlap,
                    sequence_provenance=prov,
                )
            )

        # 3p rDNA resolution
        g_chain_3p: list[SegmentPort] = []
        if last_port.name in rdna_3p:
            r = rdna_3p[last_port.name]
            has_first_edge = any(e.to_port == r.path[1] for e in graph.get_outgoing(r.path[0]))
            overlap = _find_overlap_len(graph, r.path[0], r.path[1]) if has_first_edge else 0
            if comps:
                comps[-1].overlap_to_next_bp = overlap
            if not has_first_edge:
                _add_comp(
                    PathComponent(
                        segment_id="rdna_tangle_gap",
                        strand="+",
                        is_gap=True,
                        gap_length=DEFAULT_AGP_POLICY.unknown_gap_bp.value,
                        gap_type="repeat",
                        linkage="yes",
                        linkage_evidence="rna_tangle",
                        overlap_to_next_bp=0,
                        gap_component_type="U",
                    )
                )

            for i in range(1, len(r.path)):
                port = r.path[i]
                p_seg = graph.segments[port.name]
                p_seq = p_seg.sequence if port.strand == Strand.FORWARD else reverse_complement(p_seg.sequence or "")
                nxt_port = r.path[i + 1] if i < len(r.path) - 1 else None
                has_nxt_edge = any(e.to_port == nxt_port for e in graph.get_outgoing(port)) if nxt_port else False
                nxt_ov = _find_overlap_len(graph, port, nxt_port) if has_nxt_edge else 0
                prov = (
                    SequenceProvenance.HIFI_GERMLINE
                    if (i == len(r.path) - 1 and not (r.capped_via_heuristic and r.distal_cap_port is not None))
                    else SequenceProvenance.RDNA_CHAINED
                )
                if _add_comp(
                    PathComponent(
                        segment_id=port.name,
                        strand=port.strand.to_symbol(),
                        sequence=p_seq,
                        overlap_to_next_bp=nxt_ov,
                        sequence_provenance=prov,
                    )
                ):
                    if nxt_port is not None and not has_nxt_edge:
                        _add_comp(
                            PathComponent(
                                segment_id="rdna_tangle_gap",
                                strand="+",
                                is_gap=True,
                                gap_length=DEFAULT_AGP_POLICY.unknown_gap_bp.value,
                                gap_type="repeat",
                                linkage="yes",
                                linkage_evidence="rna_tangle",
                                overlap_to_next_bp=0,
                                gap_component_type="U",
                            )
                        )

            if r.capped_via_heuristic and r.distal_cap_port is not None:
                cap_port = r.distal_cap_port
                cap_seg = graph.segments.get(cap_port.name)
                cap_strand = cap_port.strand.to_symbol()
                cap_seq = (
                    cap_seg.sequence
                    if (cap_strand == "+" and cap_seg is not None)
                    else reverse_complement(cap_seg.sequence if cap_seg is not None else "")
                )
                _add_comp(
                    PathComponent(
                        segment_id="rdna_tangle_gap",
                        strand="+",
                        is_gap=True,
                        gap_length=DEFAULT_AGP_POLICY.unknown_gap_bp.value,
                        gap_type="repeat",
                        linkage="yes",
                        linkage_evidence="rna_tangle",
                        overlap_to_next_bp=0,
                        gap_component_type="U",
                    )
                )
                _add_comp(
                    PathComponent(
                        segment_id=cap_port.name,
                        strand=cap_strand,
                        sequence=cap_seq,
                        overlap_to_next_bp=0,
                        sequence_provenance=SequenceProvenance.HIFI_GERMLINE,
                    )
                )
            consumed.update(r.excluded_unitigs)
        elif last_port.name in dropout_by_anchor:
            # Somatic dropout rescue
            dqrec, derec, dseq = dropout_by_anchor[last_port.name]
            g_sym = derec.germline_port
            if g_sym:
                g_name, g_strand = _parse_port_spec(g_sym)
                _add_comp(
                    PathComponent(
                        segment_id=f"{last_port.name}_ont_rescue",
                        strand="+",
                        sequence=dseq,
                        overlap_to_next_bp=0,
                        sequence_provenance=SequenceProvenance.ONT_DERIVED,
                        unpolished=True,
                    )
                )
                if g_name in graph.segments:
                    g_seg = graph.segments[g_name]
                    g_seq = g_seg.sequence if g_strand == Strand.FORWARD else reverse_complement(g_seg.sequence or "")
                    _add_comp(
                        PathComponent(
                            segment_id=g_name,
                            strand=g_strand.to_symbol(),
                            sequence=g_seq,
                            overlap_to_next_bp=0,
                            sequence_provenance=SequenceProvenance.HIFI_GERMLINE,
                        )
                    )
        else:
            # 3p germline continuation extension
            if last_port.to_symbol() in junctions_by_stem:
                j_3p = junctions_by_stem[last_port.to_symbol()]
                if j_3p.germline_continuation_port:
                    gp_3p = _parse_port(j_3p.germline_continuation_port)
                    g_chain_3p = _trace_germline_chain(graph, ledger_map, gp_3p, max_hops, core_coverage=core_coverage)

            if g_chain_3p:
                if comps:
                    comps[-1].overlap_to_next_bp = _find_overlap_len(graph, last_port, g_chain_3p[0])
                for i, gp in enumerate(g_chain_3p):
                    nxt_ov = _find_overlap_len(graph, gp, g_chain_3p[i + 1]) if i < len(g_chain_3p) - 1 else 0
                    gp_seg = graph.segments[gp.name]
                    gp_seq = (
                        gp_seg.sequence if gp.strand == Strand.FORWARD else reverse_complement(gp_seg.sequence or "")
                    )
                    _add_comp(
                        PathComponent(
                            segment_id=gp.name,
                            strand=gp.strand.to_symbol(),
                            sequence=gp_seq,
                            overlap_to_next_bp=nxt_ov,
                            sequence_provenance=SequenceProvenance.HIFI_GERMLINE,
                        )
                    )

        # 3p gap resolution
        cand_left = g_chain_3p[-1] if g_chain_3p else last_port
        gres = gap_by_left.get((cand_left.name, cand_left.strand.value))
        if gres is not None and gres.resolution_record and gres.resolution_record.side == "5p":
            gres = None
        if gres is None:
            cand_gres = gap_by_left_seg.get(cand_left.name)
            if cand_gres is not None and cand_gres.resolution_record and cand_gres.resolution_record.side != "5p":
                gres = cand_gres
        if gres is None and cand_left != last_port:
            cand_gres = gap_by_left.get((last_port.name, last_port.strand.value)) or gap_by_left_seg.get(last_port.name)
            if cand_gres is not None and cand_gres.resolution_record and cand_gres.resolution_record.side != "5p":
                gres = cand_gres
        if gres is None and not g_chain_3p:
            for e in graph.get_outgoing(last_port):
                target_seg = e.to_port.name
                if (target_seg, e.to_port.strand.value) in gap_by_left or target_seg in gap_by_left_seg:
                    cand_gres = gap_by_left.get((target_seg, e.to_port.strand.value)) or gap_by_left_seg.get(target_seg)
                    if (
                        cand_gres is not None
                        and cand_gres.resolution_record
                        and cand_gres.resolution_record.side != "5p"
                    ):
                        gres = cand_gres
                        e_ov = parse_dovetail_overlap_len(e.cigar) if e.cigar and e.cigar != "*" else e.overlap_len
                        if comps:
                            comps[-1].overlap_to_next_bp = e_ov
                        grs_seg = graph.segments[target_seg]
                        grs_seq = (
                            grs_seg.sequence
                            if e.to_port.strand == Strand.FORWARD
                            else reverse_complement(grs_seg.sequence or "")
                        )
                        _add_comp(
                            PathComponent(
                                segment_id=target_seg,
                                strand=e.to_port.strand.to_symbol(),
                                sequence=grs_seq,
                                overlap_to_next_bp=0,
                                sequence_provenance=SequenceProvenance.HIFI_GERMLINE,
                            )
                        )
                        break

        if gres is not None:
            qrec = gres.resolution_record
            r_name = (
                qrec.right_anchor_port.rstrip("+-")
                if (qrec and qrec.right_anchor_port)
                else (gres.right_flank_port.name if gres.right_flank_port else None)
            )
            if qrec is not None:
                placed_u: list[str] = []
                if gres.intermediate_ports:
                    for ip in gres.intermediate_ports:
                        if ip.name in graph.segments:
                            placed_u.append(ip.name)
                elif gres.representative_sequence and gres.representative_read:
                    placed_u.append(gres.representative_read.qname)
                if r_name and r_name in graph.segments:
                    placed_u.append(r_name)
                if placed_u:
                    qrec.unitig_ids = tuple(placed_u)
            r_strand = (
                Strand.FORWARD
                if (qrec and qrec.right_anchor_port and "+" in qrec.right_anchor_port)
                else (gres.right_flank_port.strand if gres.right_flank_port else Strand.FORWARD)
            )
            if gres.intermediate_ports:
                gap_ev = qrec.agp_linkage_evidence if (qrec and qrec.agp_linkage_evidence) else "satellite-kmer"
                for ip in gres.intermediate_ports:
                    _add_comp(
                        PathComponent(
                            segment_id=f"gap_{gres.chrom}_{gres.side}_{gres.gap_index}_anchor",
                            strand="+",
                            is_gap=True,
                            gap_length=DEFAULT_AGP_POLICY.unknown_gap_bp.value,
                            gap_type="repeat",
                            linkage="yes",
                            linkage_evidence=gap_ev,
                            overlap_to_next_bp=0,
                            gap_component_type="U",
                        )
                    )
                    if ip.name in graph.segments:
                        ip_seg = graph.segments[ip.name]
                        ip_seq = (
                            ip_seg.sequence
                            if ip.strand == Strand.FORWARD
                            else reverse_complement(ip_seg.sequence or "")
                        )
                        _add_comp(
                            PathComponent(
                                segment_id=ip.name,
                                strand=ip.strand.to_symbol(),
                                sequence=ip_seq,
                                overlap_to_next_bp=0,
                                sequence_provenance=SequenceProvenance.HIFI_GERMLINE,
                            )
                        )
                _add_comp(
                    PathComponent(
                        segment_id=f"gap_{gres.chrom}_{gres.side}_{gres.gap_index}_cap",
                        strand="+",
                        is_gap=True,
                        gap_length=DEFAULT_AGP_POLICY.unknown_gap_bp.value,
                        gap_type="repeat",
                        linkage="yes",
                        linkage_evidence=gap_ev,
                        overlap_to_next_bp=0,
                        gap_component_type="U",
                    )
                )
            elif gres.representative_sequence and gres.representative_read:
                _add_comp(
                    PathComponent(
                        segment_id=gres.representative_read.qname,
                        strand="+",
                        sequence=gres.representative_sequence,
                        overlap_to_next_bp=0,
                        sequence_provenance=SequenceProvenance.ONT_DERIVED,
                        unpolished=True,
                    )
                )
            else:
                gap_len = (
                    int(round(gres.median_gap_bp)) if gres.median_gap_bp else DEFAULT_AGP_POLICY.unknown_gap_bp.value
                )
                gap_type = (
                    "repeat"
                    if (
                        qrec
                        and qrec.resolution_type
                        in (ResolutionType.SAT_GAP_VARIABLE, ResolutionType.SUBTELOMERIC_SATELLITE_CAP)
                    )
                    else "scaffold"
                )
                gap_ev = (
                    qrec.agp_linkage_evidence
                    if (qrec and qrec.agp_linkage_evidence)
                    else ("satellite-kmer" if gap_type == "repeat" else "hifi-spanning")
                )
                comp_type = "U" if gap_len == DEFAULT_AGP_POLICY.unknown_gap_bp.value else "N"
                _add_comp(
                    PathComponent(
                        segment_id=f"gap_{gres.chrom}_{gres.side}_{gres.gap_index}",
                        strand="+",
                        is_gap=True,
                        gap_length=gap_len,
                        gap_type=gap_type,
                        linkage="yes",
                        linkage_evidence=gap_ev,
                        overlap_to_next_bp=0,
                        gap_component_type=comp_type,
                    )
                )
            if r_name and r_name in graph.segments:
                r_seg = graph.segments[r_name]
                r_seq = r_seg.sequence if r_strand == Strand.FORWARD else reverse_complement(r_seg.sequence or "")
                _add_comp(
                    PathComponent(
                        segment_id=r_name,
                        strand=r_strand.to_symbol(),
                        sequence=r_seq,
                        overlap_to_next_bp=0,
                        sequence_provenance=SequenceProvenance.HIFI_GERMLINE,
                    )
                )

        if not comps:
            continue

        scaff_id = f"{prefix}{seg_name}" if prefix else seg_name
        seq, agp_recs = assemble_scaffold_sequence(comps, scaffold_id=scaff_id)
        scaffolds[scaff_id] = seq
        all_agp.extend(agp_recs)

    if active_resolutions is not None:
        component_to_scaf: dict[str, str] = {}
        component_spans: dict[tuple[str, str], list[tuple[int, int]]] = defaultdict(list)
        for r_agp in all_agp:
            scaf_id = r_agp.object_id
            beg = r_agp.object_beg
            end = r_agp.object_end
            if r_agp.component_type == "W":
                c_id = r_agp.component_id
                component_to_scaf[c_id] = scaf_id
                component_spans[(scaf_id, c_id)].append((beg, end))

        retained: list[ResolutionRecord] = []
        for rec in list(active_resolutions):
            if (
                rec.resolution_type == ResolutionType.SAT_GAP_VARIABLE
                and rec.agp_linkage_evidence == "unbridged"
                and (rec.read_count is None or rec.read_count == 0)
                and rec.spanning_read_count is None
            ):
                continue
            target_scaf: str | None = None
            if rec.left_anchor_port:
                l_u = rec.left_anchor_port.rstrip("+-")
                if l_u in component_to_scaf:
                    target_scaf = component_to_scaf[l_u]
            if target_scaf is None and rec.chrom in scaffolds:
                target_scaf = rec.chrom

            spans: list[tuple[int, int]] = []
            placed_uids: list[str] = []

            if rec.resolution_type == ResolutionType.PDE_TERMINAL_CAP and target_scaf is not None:
                scaf_comps = [r for r in all_agp if r.object_id == target_scaf and r.component_type == "W"]
                stem_u = rec.left_anchor_port.rstrip("+-") if rec.left_anchor_port else rec.chrom
                comp_ids = [c.component_id for c in scaf_comps]
                if stem_u in comp_ids:
                    stem_idx = comp_ids.index(stem_u)
                    if rec.side == "5p" and stem_idx > 0:
                        distal_comps = scaf_comps[:stem_idx]
                        if any(c.component_id in rec.unitig_ids for c in distal_comps):
                            distal_ids = {c.component_id for c in distal_comps}
                            # Preserve chain order from rec.unitig_ids where available, prepend any new distal unitigs
                            new_distal = [c.component_id for c in distal_comps if c.component_id not in rec.unitig_ids]
                            existing_kept = [u for u in rec.unitig_ids if u in distal_ids]
                            placed_uids = new_distal + existing_kept
                            spans = [(c.object_beg, c.object_end) for c in distal_comps]
                            rec.unitig_ids = tuple(placed_uids)
                            if (
                                new_distal
                                or not rec.right_anchor_port
                                or rec.right_anchor_port.rstrip("+-") not in distal_ids
                            ):
                                rec.right_anchor_port = f"{distal_comps[0].component_id}-"
                    elif rec.side == "3p" and stem_idx < len(scaf_comps) - 1:
                        distal_comps = scaf_comps[stem_idx + 1 :]
                        if any(c.component_id in rec.unitig_ids for c in distal_comps):
                            distal_ids = {c.component_id for c in distal_comps}
                            existing_kept = [u for u in rec.unitig_ids if u in distal_ids]
                            new_distal = [c.component_id for c in distal_comps if c.component_id not in rec.unitig_ids]
                            placed_uids = existing_kept + new_distal
                            spans = [(c.object_beg, c.object_end) for c in distal_comps]
                            rec.unitig_ids = tuple(placed_uids)
                            if (
                                new_distal
                                or not rec.right_anchor_port
                                or rec.right_anchor_port.rstrip("+-") not in distal_ids
                            ):
                                rec.right_anchor_port = f"{distal_comps[-1].component_id}{distal_comps[-1].orientation}"

            if not spans:
                for u in rec.unitig_ids:
                    if (target_scaf, u) in component_spans:
                        spans.extend(component_spans[(target_scaf, u)])
                        placed_uids.append(u)
                    elif target_scaf is None and u in component_to_scaf:
                        s_id = component_to_scaf[u]
                        target_scaf = s_id
                        spans.extend(component_spans[(s_id, u)])
                        placed_uids.append(u)

            if target_scaf is None and rec.right_anchor_port:
                r_u = rec.right_anchor_port.rstrip("+-")
                if r_u in component_to_scaf:
                    target_scaf = component_to_scaf[r_u]

            if target_scaf is not None and spans:
                rec.chrom = target_scaf
                rec.unitig_ids = tuple(placed_uids)
                rec.agp_start = min(s[0] for s in spans)
                rec.agp_end = max(s[1] for s in spans)
                scaf_len = len(scaffolds[target_scaf]) if target_scaf in scaffolds else None
                if scaf_len is not None and rec.side in ("5p", "3p"):
                    if rec.agp_end == scaf_len and rec.agp_start > 1:
                        rec.side = "3p"
                    elif rec.agp_start == 1 and rec.agp_end < scaf_len:
                        rec.side = "5p"
                    rec.feature_id = build_feature_id(
                        rec.resolution_type,
                        rec.left_anchor_port.rstrip("+-") if rec.left_anchor_port else rec.chrom,
                        rec.side,
                        index=1,
                    )
                if rec.resolution_type == ResolutionType.INTERCALARY_GRS_BRIDGE:
                    if rec.left_anchor_port and rec.right_anchor_port:
                        l_u = rec.left_anchor_port.rstrip("+-")
                        r_u = rec.right_anchor_port.rstrip("+-")
                        l_spans = component_spans.get((target_scaf, l_u))
                        r_spans = component_spans.get((target_scaf, r_u))
                        if l_spans and r_spans:
                            if min(s[0] for s in l_spans) > min(s[0] for s in r_spans):
                                rec.left_anchor_port, rec.right_anchor_port = (
                                    rec.right_anchor_port,
                                    rec.left_anchor_port,
                                )
                                rec.excluded_unitig_ids = tuple(reversed(rec.excluded_unitig_ids))
                                rec.feature_id = build_feature_id(
                                    rec.resolution_type, rec.left_anchor_port.rstrip("+-"), rec.side, index=1
                                )
                retained.append(rec)
            elif target_scaf is not None and not spans and rec.resolution_type != ResolutionType.SAT_GAP_VARIABLE:
                other_scafs = {
                    u: component_to_scaf[u]
                    for u in rec.unitig_ids
                    if u in component_to_scaf and component_to_scaf[u] != target_scaf and component_to_scaf[u] != u
                }
                filtered_uids = tuple(u for u in rec.unitig_ids if u not in other_scafs)
                if filtered_uids and filtered_uids != rec.unitig_ids:
                    rec.unitig_ids = filtered_uids
                    if rec.right_anchor_port and rec.right_anchor_port.rstrip("+-") not in filtered_uids:
                        rec.right_anchor_port = f"{filtered_uids[-1]}-"
                anchor_spans: list[tuple[int, int]] = []
                if rec.left_anchor_port:
                    l_u = rec.left_anchor_port.rstrip("+-")
                    if (target_scaf, l_u) in component_spans:
                        anchor_spans.extend(component_spans[(target_scaf, l_u)])
                if rec.right_anchor_port:
                    r_u = rec.right_anchor_port.rstrip("+-")
                    if (target_scaf, r_u) in component_spans:
                        anchor_spans.extend(component_spans[(target_scaf, r_u)])
                if anchor_spans:
                    rec.chrom = target_scaf
                    rec.agp_start = min(s[0] for s in anchor_spans)
                    rec.agp_end = max(s[1] for s in anchor_spans)
                    scaf_len = len(scaffolds[target_scaf]) if target_scaf in scaffolds else None
                    if scaf_len is not None and rec.side in ("5p", "3p"):
                        if rec.agp_end == scaf_len and rec.agp_start > 1:
                            rec.side = "3p"
                        elif rec.agp_start == 1:
                            rec.side = "5p"
                        rec.feature_id = build_feature_id(
                            rec.resolution_type,
                            rec.left_anchor_port.rstrip("+-") if rec.left_anchor_port else rec.chrom,
                            rec.side,
                            index=1,
                        )
                    retained.append(rec)

        active_resolutions.clear()
        active_resolutions.extend(retained)

    if nominate_super:
        scaffolds, all_agp, _ = nominate_super_scaffolds(
            scaffolds=scaffolds,
            agp_records=all_agp,
            junctions=junctions,
            min_super_length=min_super_length,
            resolutions=active_resolutions,
            end_states=list(end_states) if isinstance(end_states, list) else None,
        )

    return scaffolds, all_agp


def run_report(args: argparse.Namespace) -> None:
    """Execute the report subcommand."""
    graph = load_gfa(args.gfa)
    annotations = aggregate_node_annotations(
        tiara=args.tiara,
        barrnap=args.barrnap,
        sequence_ids=graph.segments.keys(),
    )
    ledger = build_end_evidence_ledger(
        graph=graph,
        annotations=annotations,
        terminal_window_bp=args.terminal_window,
    )
    if args.out is not None:
        emit_tsv_ledger(ledger, args.out)
    else:
        emit_tsv_ledger(ledger, sys.stdout)


def run_finish(args: argparse.Namespace) -> None:
    """Execute the finish subcommand."""
    args.out_dir.mkdir(parents=True, exist_ok=True)
    graph = load_gfa(args.gfa, load_sequences=True)
    annotations = aggregate_node_annotations(
        tiara=args.tiara,
        barrnap=args.barrnap,
        sequence_ids=graph.segments.keys(),
    )
    ledger = build_end_evidence_ledger(
        graph=graph,
        annotations=annotations,
        terminal_window_bp=args.terminal_window,
    )
    ledger_map = {ev.port.to_symbol(): ev for ev in ledger} | {f"{ev.segment_name}:{ev.end}": ev for ev in ledger}
    core_coverage = compute_core_coverage(graph, ledger)
    junction_records, end_state_records, resolution_records = collect_features(
        graph=graph,
        ledger=ledger,
        core_coverage=core_coverage,
    )

    engine = OntBridgingEngine(graph, ledger_map)
    orphan_caps = engine.identify_orphan_caps()

    if args.ont_paf is not None:
        pre_somatic: dict[tuple[str, str], tuple[str | None, float | None]] = {
            (rec.chrom, rec.side): (rec.somatic_port, rec.somatic_depth_ratio) for rec in end_state_records
        }

        end_state_records, resolution_records = engine.identify_and_resolve_gaps(
            paf_source=args.ont_paf,
            end_states=end_state_records,
            quirks=resolution_records,
            reads_sequences=args.ont_reads,
        )

        for rec in end_state_records:
            if (
                rec.state == EndState.END_SOMATIC_MISSING
                and rec.somatic_port
                and rec.germline_port
                and not engine._has_terminal_telomere(rec.germline_port)
            ):
                stem_p = SegmentPort(rec.chrom, Strand.FORWARD if rec.side == "3p" else Strand.REVERSE)
                g_name, g_strand = _parse_port_spec(rec.germline_port)
                s_name, s_strand = _parse_port_spec(rec.somatic_port)
                d_qrec, d_erec, d_seq = engine.resolve_somatic_dropout(
                    chrom=rec.chrom,
                    side=rec.side,
                    stem_port=stem_p,
                    germline_port=SegmentPort(g_name, g_strand),
                    somatic_port=SegmentPort(s_name, s_strand),
                    paf_path=args.ont_paf,
                    reads_fasta=args.ont_reads,
                    core_coverage_bp_median=rec.core_coverage_bp_median,
                    somatic_depth_ratio=rec.somatic_depth_ratio,
                    germline_depth_ratio=rec.germline_depth_ratio,
                )
                if d_erec.resolved:
                    rec.state = d_erec.state
                    rec.resolved = d_erec.resolved
                    rec.notes = d_erec.notes
                    resolution_records.append(d_qrec)

        for rec in end_state_records:
            som_port, som_ratio = pre_somatic.get((rec.chrom, rec.side), (None, None))
            if som_port is not None and rec.somatic_port != som_port:
                rec.somatic_port = som_port
            if som_ratio is not None and rec.somatic_depth_ratio is None:
                rec.somatic_depth_ratio = som_ratio
            if rec.germline_port is not None:
                cap_seg = rec.germline_port.rstrip("+-")
                if cap_seg in orphan_caps:
                    rec.germline_port = orphan_caps[cap_seg][0]

    orphan_names = set(orphan_caps.keys())
    end_state_records = [rec for rec in end_state_records if rec.chrom not in orphan_names]

    # rDNA cassette resolution
    rdna_resolver = RdnaCassetteResolver(graph, ledger_map, annotations=annotations)
    entry_ports = rdna_resolver.find_rdna_entry_ports()
    rdna_resolutions: list[RdnaResolution] = []
    rdna_chained_names: set[str] = set()
    consumed_rdna: set[str] = set()

    for ep in entry_ports:
        if ep.name in rdna_chained_names or ep.name in consumed_rdna:
            continue
        chrom = ep.name
        side = "3p" if ep.strand == Strand.FORWARD else "5p"
        res = rdna_resolver.resolve_rdna_cassette(chrom=chrom, side=side, entry_port=ep, consumed_rdna=consumed_rdna)
        if len(res.path) > 1:
            rdna_resolutions.append(res)
            for p in res.path[1:]:
                rdna_chained_names.add(p.name)
                consumed_rdna.add(p.name)
            consumed_rdna.update(res.excluded_unitigs)
            if res.distal_cap_port is not None:
                rdna_chained_names.add(res.distal_cap_port.name)
                consumed_rdna.add(res.distal_cap_port.name)

            if res.reached_distal_cap or res.capped_via_heuristic or res.is_intercalary:
                rec_to_add = res.end_state_record
                if res.capped_via_heuristic and res.distal_cap_port is not None:
                    rec_to_add.germline_port = res.distal_cap_port.to_symbol()
                    rec_to_add.resolved = True
                elif res.is_intercalary:
                    rec_to_add.resolved = True

                new_ends: list[EndStateRecord] = []
                replaced = False
                for er in end_state_records:
                    if er.chrom == chrom and er.side == side:
                        new_ends.append(rec_to_add)
                        replaced = True
                    else:
                        new_ends.append(er)
                if not replaced:
                    new_ends.append(rec_to_add)
                end_state_records = new_ends

            new_res: list[ResolutionRecord] = []
            rep_r = False
            for qr in resolution_records:
                if qr.chrom == chrom and qr.side == side and qr.resolution_type == ResolutionType.RDNA_CASSETTE:
                    new_res.append(res.resolution_record)
                    rep_r = True
                else:
                    new_res.append(qr)
            if not rep_r:
                new_res.append(res.resolution_record)
            resolution_records = new_res

    end_state_records = [rec for rec in end_state_records if rec.chrom not in rdna_chained_names]

    # Assemble scaffolds and AGP records
    prefix = getattr(args, "prefix", "") or ""
    nominate_super = getattr(args, "nominate_super", False)
    min_super_length = getattr(args, "min_super_length", DEFAULT_MIN_SUPER_LENGTH)
    scaffolds, agp_records = build_scaffolds(
        graph=graph,
        ledger_map=ledger_map,
        end_states=end_state_records,
        resolutions=resolution_records,
        gap_resolutions=engine.gap_resolutions,
        dropout_resolutions=engine.dropout_resolutions,
        rdna_resolutions=rdna_resolutions,
        prefix=prefix,
        junctions=junction_records,
        nominate_super=nominate_super,
        min_super_length=min_super_length,
        core_coverage=core_coverage,
    )

    fasta_target: Path
    if args.emit_fasta is not None and not isinstance(args.emit_fasta, bool):
        fasta_target = Path(args.emit_fasta)
    else:
        fasta_target = args.out_dir / "scaffolds.fasta"

    agp_target: Path
    if args.emit_agp is not None and not isinstance(args.emit_agp, bool):
        agp_target = Path(args.emit_agp)
    else:
        agp_target = args.out_dir / "scaffolds.agp"

    emit_fasta(scaffolds, fasta_target)
    emit_agp(agp_records, agp_target)

    # Emit eliminated GRS regions BED
    grs_target: Path
    if getattr(args, "emit_grs_bed", None) is not None and not isinstance(args.emit_grs_bed, bool):
        grs_target = Path(args.emit_grs_bed)
    else:
        grs_target = args.out_dir / "GRS.bed"
    emit_grs_bed(agp_records, grs_target)

    # Emit germline terminal telomeres BED
    telomere_target: Path
    if getattr(args, "emit_telomere_bed", None) is not None and not isinstance(args.emit_telomere_bed, bool):
        telomere_target = Path(args.emit_telomere_bed)
    else:
        telomere_target = args.out_dir / "telomeres.bed"
    emit_telomere_bed(scaffolds, telomere_target, policy=DEFAULT_POLICY)

    # Compute and emit summary metrics on break-site verified sequences
    summary_target: Path
    if getattr(args, "summary", None) is not None and not isinstance(args.summary, bool):
        summary_target = Path(args.summary)
    else:
        summary_target = args.out_dir / "summary.tsv"
    summary_data = compute_break_site_summary(
        graph=graph,
        scaffolds=scaffolds,
        agp_records=agp_records,
        junctions=junction_records,
        resolutions=resolution_records,
    )
    emit_summary_tsv(summary_data, summary_target)
    print(format_summary_table(summary_data))

    excluded_target: Path
    if args.excluded_fasta is not None and not isinstance(args.excluded_fasta, bool):
        excluded_target = Path(args.excluded_fasta)
    else:
        excluded_target = args.out_dir / "excluded.fasta"

    placed_uids = {str(r.component_id) for r in agp_records if getattr(r, "component_type", "") == "W"}
    emit_excluded_fasta(graph, resolution_records, excluded_target, placed_unitigs=placed_uids)

    emit_tsv_ledger(ledger, args.out_dir / "ledger.tsv")
    emit_junctions_tsv(junction_records, args.out_dir / "junctions.tsv")
    emit_end_states_tsv(end_state_records, args.out_dir / "end_states.tsv")
    res_target = (
        args.resolutions if getattr(args, "resolutions", None) is not None else args.out_dir / "resolutions.tsv"
    )
    emit_resolutions_tsv(resolution_records, res_target)


def main(args: list[str] | None = None) -> int:
    """Entrypoint for the gfapdefinisher CLI."""
    parser = create_parser()
    parsed = parser.parse_args(args)

    if parsed.command is None:
        parser.print_help()
        sys.exit(0)

    if parsed.command == "report":
        run_report(parsed)
        return 0

    if parsed.command == "finish":
        run_finish(parsed)
        return 0

    return 0


if __name__ == "__main__":
    sys.exit(main())
