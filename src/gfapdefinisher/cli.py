"""Command-line interface for gfapdefinisher."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from gfapdefinisher.annotations import aggregate_node_annotations
from gfapdefinisher.classifier import compute_core_coverage
from gfapdefinisher.graph import load_gfa
from gfapdefinisher.ledger import build_end_evidence_ledger, emit_tsv_ledger
from gfapdefinisher.quirks import (
    collect_features,
    emit_end_states_tsv,
    emit_junctions_tsv,
    emit_quirks_tsv,
)
from gfapdefinisher.resolver import OntBridgingEngine, load_reads_sequences


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

    return parser


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
    core_coverage = compute_core_coverage(graph, ledger)
    junction_records, end_state_records, quirk_records = collect_features(
        graph=graph,
        ledger=ledger,
        core_coverage=core_coverage,
    )

    engine = OntBridgingEngine(graph, ledger)
    orphan_caps = engine.identify_orphan_caps()

    if args.ont_paf is not None:
        reads_dict = load_reads_sequences(args.ont_reads) if args.ont_reads is not None else None
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
