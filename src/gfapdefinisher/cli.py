"""Command-line interface for gfapdefinisher."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from gfapdefinisher.annotations import aggregate_node_annotations
from gfapdefinisher.graph import load_gfa
from gfapdefinisher.ledger import build_end_evidence_ledger, emit_tsv_ledger


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

    return 0


if __name__ == "__main__":
    sys.exit(main())
