"""Command-line interface for gfapdefinisher."""

from __future__ import annotations

import argparse
import sys
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from gfapdefinisher.annotations import aggregate_node_annotations
from gfapdefinisher.cigar import cigar_query_consumption
from gfapdefinisher.classifier import compute_core_coverage
from gfapdefinisher.emitter import (
    AgpRecord,
    PathComponent,
    SequenceProvenance,
    assemble_scaffold_sequence,
    emit_agp,
    emit_fasta,
)
from gfapdefinisher.graph import BidirectedGraph, load_gfa
from gfapdefinisher.ledger import build_end_evidence_ledger, emit_tsv_ledger
from gfapdefinisher.models import SegmentPort, Strand
from gfapdefinisher.quirks import (
    EndState,
    EndStateRecord,
    IntervalQuirk,
    QuirkRecord,
    collect_features,
    emit_end_states_tsv,
    emit_junctions_tsv,
    emit_quirks_tsv,
)
from gfapdefinisher.rdna import RdnaCassetteResolver, RdnaResolution
from gfapdefinisher.resolver import (
    OntBridgingEngine,
    SatelliteGapResolution,
    _parse_port_spec,
    load_reads_sequences,
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
    report_parser.add_argument("--tidk", type=Path, default=None, help="Path to TIDK telomeric repeat windows TSV")
    report_parser.add_argument("--barrnap", type=Path, default=None, help="Path to Barrnap rDNA annotations GFF3")
    report_parser.add_argument("--compleasm", type=Path, default=None, help="Path to Compleasm full table TSV")
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
    finish_parser.add_argument("--tidk", type=Path, default=None, help="Path to TIDK telomeric repeat windows TSV")
    finish_parser.add_argument("--barrnap", type=Path, default=None, help="Path to Barrnap rDNA annotations GFF3")
    finish_parser.add_argument("--compleasm", type=Path, default=None, help="Path to Compleasm full table TSV")
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
        "--prefix",
        type=str,
        default="",
        help="Optional prefix for scaffold identifiers",
    )

    return parser


def build_scaffolds(
    graph: BidirectedGraph,
    ledger_map: Mapping[str, Any],
    end_states: Sequence[EndStateRecord],
    quirks: Sequence[QuirkRecord],
    gap_resolutions: Sequence[SatelliteGapResolution],
    dropout_resolutions: Sequence[tuple[QuirkRecord, EndStateRecord, str | None]],
    rdna_resolutions: Sequence[RdnaResolution],
    prefix: str = "",
) -> tuple[dict[str, str], list[AgpRecord]]:
    """Assemble scaffold sequences and AGP v2.1 records with empirical provenance."""
    rdna_5p: dict[str, RdnaResolution] = {}
    rdna_3p: dict[str, RdnaResolution] = {}
    for r in rdna_resolutions:
        if len(r.path) > 1:
            if r.side == "5p":
                rdna_5p[r.path[0].name] = r
            else:
                rdna_3p[r.path[0].name] = r

    gap_by_left: dict[str, SatelliteGapResolution] = {}
    for g in gap_resolutions:
        if g.representative_sequence and g.representative_read:
            q = g.quirk_record
            left_name = q.left_anchor_port.rstrip("+-") if q.left_anchor_port else g.chrom
            gap_by_left[left_name] = g

    dropout_by_anchor: dict[str, tuple[QuirkRecord, EndStateRecord, str]] = {}
    for qrec, erec, dseq in dropout_resolutions:
        if dseq:
            dropout_by_anchor[erec.chrom] = (qrec, erec, dseq)

    consumed: set[str] = set()
    for r in rdna_resolutions:
        if len(r.path) > 1:
            for p in r.path[1:]:
                consumed.add(p.name)

    scaffolds: dict[str, str] = {}
    all_agp: list[AgpRecord] = []

    sorted_segs = sorted(graph.segments.keys(), key=lambda s: (-graph.segments[s].length, s))

    for seg_name in sorted_segs:
        if seg_name in consumed:
            continue
        seg = graph.segments[seg_name]
        if not seg.sequence:
            continue

        comps: list[PathComponent] = []

        # 5p rDNA resolution
        if seg_name in rdna_5p:
            r = rdna_5p[seg_name]
            twin_path = [p.invert() for p in reversed(r.path)]
            for i, port in enumerate(twin_path):
                p_seg = graph.segments[port.name]
                p_seq = p_seg.sequence if port.strand == Strand.FORWARD else reverse_complement(p_seg.sequence or "")
                overlap = 0
                if i < len(twin_path) - 1:
                    nxt = twin_path[i + 1]
                    for e in graph.get_outgoing(port):
                        if e.to_port == nxt:
                            overlap = cigar_query_consumption(e.cigar) if e.cigar and e.cigar != "*" else e.overlap_len
                            break
                prov = (
                    SequenceProvenance.RDNA_CHAINED
                    if (
                        "rdna" in port.name.lower() or "rrna" in port.name.lower() or (i > 0 and i < len(twin_path) - 1)
                    )
                    else SequenceProvenance.HIFI_CORE
                )
                if i == 0:
                    prov = SequenceProvenance.HIFI_GERMLINE
                comps.append(
                    PathComponent(
                        segment_id=port.name,
                        strand=port.strand.to_symbol(),
                        sequence=p_seq,
                        overlap_to_next_bp=overlap,
                        sequence_provenance=prov,
                    )
                )
                if port.name != seg_name:
                    consumed.add(port.name)
        else:
            comps.append(
                PathComponent(
                    segment_id=seg_name,
                    strand="+",
                    sequence=seg.sequence,
                    overlap_to_next_bp=0,
                    sequence_provenance=SequenceProvenance.HIFI_CORE,
                )
            )

        # 3p rDNA resolution
        if seg_name in rdna_3p:
            r = rdna_3p[seg_name]
            overlap = 0
            for e in graph.get_outgoing(r.path[0]):
                if e.to_port == r.path[1]:
                    overlap = cigar_query_consumption(e.cigar) if e.cigar and e.cigar != "*" else e.overlap_len
                    break
            comps[-1].overlap_to_next_bp = overlap

            for i in range(1, len(r.path)):
                port = r.path[i]
                p_seg = graph.segments[port.name]
                p_seq = p_seg.sequence if port.strand == Strand.FORWARD else reverse_complement(p_seg.sequence or "")
                nxt_ov = 0
                if i < len(r.path) - 1:
                    nxt = r.path[i + 1]
                    for e in graph.get_outgoing(port):
                        if e.to_port == nxt:
                            nxt_ov = cigar_query_consumption(e.cigar) if e.cigar and e.cigar != "*" else e.overlap_len
                            break
                prov = SequenceProvenance.HIFI_GERMLINE if i == len(r.path) - 1 else SequenceProvenance.RDNA_CHAINED
                comps.append(
                    PathComponent(
                        segment_id=port.name,
                        strand=port.strand.to_symbol(),
                        sequence=p_seq,
                        overlap_to_next_bp=nxt_ov,
                        sequence_provenance=prov,
                    )
                )
                consumed.add(port.name)

        # 3p gap resolution
        gres = gap_by_left.get(seg_name)
        if gres is None:
            for e in graph.get_outgoing(SegmentPort(seg_name, Strand.FORWARD)):
                target_seg = e.to_port.name
                if target_seg in gap_by_left:
                    cand_gres = gap_by_left[target_seg]
                    gres = cand_gres
                    e_ov = cigar_query_consumption(e.cigar) if e.cigar and e.cigar != "*" else e.overlap_len
                    comps[-1].overlap_to_next_bp = e_ov
                    grs_seg = graph.segments[target_seg]
                    grs_seq = (
                        grs_seg.sequence
                        if e.to_port.strand == Strand.FORWARD
                        else reverse_complement(grs_seg.sequence or "")
                    )
                    comps.append(
                        PathComponent(
                            segment_id=target_seg,
                            strand=e.to_port.strand.to_symbol(),
                            sequence=grs_seq,
                            overlap_to_next_bp=0,
                            sequence_provenance=SequenceProvenance.HIFI_GERMLINE,
                        )
                    )
                    consumed.add(target_seg)
                    break

        if gres is not None and gres.representative_sequence and gres.representative_read:
            qrec = gres.quirk_record
            r_name = qrec.right_anchor_port.rstrip("+-") if qrec.right_anchor_port else None
            r_strand = Strand.FORWARD if (qrec.right_anchor_port and "+" in qrec.right_anchor_port) else Strand.REVERSE
            comps.append(
                PathComponent(
                    segment_id=gres.representative_read.qname,
                    strand="+",
                    sequence=gres.representative_sequence,
                    overlap_to_next_bp=0,
                    sequence_provenance=SequenceProvenance.ONT_DERIVED,
                    unpolished=True,
                )
            )
            if r_name and r_name in graph.segments:
                r_seg = graph.segments[r_name]
                r_seq = r_seg.sequence if r_strand == Strand.FORWARD else reverse_complement(r_seg.sequence or "")
                comps.append(
                    PathComponent(
                        segment_id=r_name,
                        strand=r_strand.to_symbol(),
                        sequence=r_seq,
                        overlap_to_next_bp=0,
                        sequence_provenance=SequenceProvenance.HIFI_GERMLINE,
                    )
                )
                consumed.add(r_name)

        # Somatic dropout rescue
        if seg_name in dropout_by_anchor:
            dqrec, derec, dseq = dropout_by_anchor[seg_name]
            g_sym = derec.germline_port
            if g_sym:
                g_name, g_strand = _parse_port_spec(g_sym)
                comps.append(
                    PathComponent(
                        segment_id=f"{seg_name}_ont_rescue",
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
                    comps.append(
                        PathComponent(
                            segment_id=g_name,
                            strand=g_strand.to_symbol(),
                            sequence=g_seq,
                            overlap_to_next_bp=0,
                            sequence_provenance=SequenceProvenance.HIFI_GERMLINE,
                        )
                    )
                    consumed.add(g_name)

        scaff_id = f"{prefix}{seg_name}" if prefix else seg_name
        seq, agp_recs = assemble_scaffold_sequence(comps, scaffold_id=scaff_id)
        scaffolds[scaff_id] = seq
        all_agp.extend(agp_recs)
        consumed.add(seg_name)

    return scaffolds, all_agp


def run_report(args: argparse.Namespace) -> None:
    """Execute the report subcommand."""
    graph = load_gfa(args.gfa)
    annotations = aggregate_node_annotations(
        tiara=args.tiara,
        tidk=args.tidk,
        barrnap=args.barrnap,
        compleasm=args.compleasm,
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
        tidk=args.tidk,
        barrnap=args.barrnap,
        compleasm=args.compleasm,
        sequence_ids=graph.segments.keys(),
    )
    ledger = build_end_evidence_ledger(
        graph=graph,
        annotations=annotations,
        terminal_window_bp=args.terminal_window,
    )
    ledger_map = {ev.port.to_symbol(): ev for ev in ledger} | {f"{ev.segment_name}:{ev.end}": ev for ev in ledger}
    core_coverage = compute_core_coverage(graph, ledger)
    junction_records, end_state_records, quirk_records = collect_features(
        graph=graph,
        ledger=ledger,
        core_coverage=core_coverage,
    )

    engine = OntBridgingEngine(graph, ledger_map)
    orphan_caps = engine.identify_orphan_caps()

    reads_dict: dict[str, str] | None = None
    if args.ont_paf is not None:
        if args.ont_reads is not None:
            target_ids = engine.collect_spanning_read_ids(args.ont_paf)
            for er in end_state_records:
                if er.state == EndState.END_SOMATIC_MISSING and er.somatic_port and er.germline_port:
                    stem_p = SegmentPort(er.chrom, Strand.FORWARD if er.side == "3p" else Strand.REVERSE)
                    g_name, g_strand = _parse_port_spec(er.germline_port)
                    for sr in engine.collect_spanning_reads(args.ont_paf, stem_p, SegmentPort(g_name, g_strand)):
                        target_ids.add(sr.qname)
            reads_dict = load_reads_sequences(args.ont_reads, target_qnames=target_ids)

        pre_somatic: dict[tuple[str, str], tuple[str | None, float | None]] = {
            (rec.chrom, rec.side): (rec.somatic_port, rec.somatic_depth_ratio) for rec in end_state_records
        }

        end_state_records, quirk_records = engine.identify_and_resolve_gaps(
            paf_source=args.ont_paf,
            end_states=end_state_records,
            quirks=quirk_records,
            reads_sequences=reads_dict,
        )

        for rec in end_state_records:
            if rec.state == EndState.END_SOMATIC_MISSING and rec.somatic_port and rec.germline_port:
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
                    reads_fasta=reads_dict,
                    core_coverage_bp_median=rec.core_coverage_bp_median,
                    somatic_depth_ratio=rec.somatic_depth_ratio,
                    germline_depth_ratio=rec.germline_depth_ratio,
                )
                if d_erec.resolved:
                    rec.state = d_erec.state
                    rec.resolved = d_erec.resolved
                    rec.notes = d_erec.notes
                    quirk_records.append(d_qrec)

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

    for ep in entry_ports:
        if ep.name in rdna_chained_names:
            continue
        chrom = ep.name
        side = "3p" if ep.strand == Strand.FORWARD else "5p"
        res = rdna_resolver.resolve_rdna_cassette(chrom=chrom, side=side, entry_port=ep)
        if len(res.path) > 1:
            rdna_resolutions.append(res)
            for p in res.path[1:]:
                rdna_chained_names.add(p.name)

            new_ends: list[EndStateRecord] = []
            replaced = False
            for er in end_state_records:
                if er.chrom == chrom and er.side == side:
                    new_ends.append(res.end_state_record)
                    replaced = True
                else:
                    new_ends.append(er)
            if not replaced:
                new_ends.append(res.end_state_record)
            end_state_records = new_ends

            new_quirks: list[QuirkRecord] = []
            rep_q = False
            for qr in quirk_records:
                if qr.chrom == chrom and qr.side == side and qr.quirk_type == IntervalQuirk.RDNA_CASSETTE:
                    new_quirks.append(res.quirk_record)
                    rep_q = True
                else:
                    new_quirks.append(qr)
            if not rep_q:
                new_quirks.append(res.quirk_record)
            quirk_records = new_quirks

    end_state_records = [rec for rec in end_state_records if rec.chrom not in rdna_chained_names]

    # Assemble scaffolds and AGP records
    prefix = getattr(args, "prefix", "") or ""
    scaffolds, agp_records = build_scaffolds(
        graph=graph,
        ledger_map=ledger_map,
        end_states=end_state_records,
        quirks=quirk_records,
        gap_resolutions=engine.gap_resolutions,
        dropout_resolutions=engine.dropout_resolutions,
        rdna_resolutions=rdna_resolutions,
        prefix=prefix,
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

    emit_tsv_ledger(ledger, args.out_dir / "ledger.tsv")
    emit_junctions_tsv(junction_records, args.out_dir / "junctions.tsv")
    emit_end_states_tsv(end_state_records, args.out_dir / "end_states.tsv")
    emit_quirks_tsv(quirk_records, args.out_dir / "quirks.tsv")


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
