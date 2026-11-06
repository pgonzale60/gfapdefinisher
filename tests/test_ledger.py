"""Unit and integration tests for unitig end-evidence ledger generation."""

import io
from dataclasses import FrozenInstanceError
from pathlib import Path

import pytest

from gfapdefinisher.annotations import (
    BarrnapHit,
    CompleasmHit,
    NodeAnnotations,
    TiaraHit,
    TidkWindow,
)
from gfapdefinisher.cli import main
from gfapdefinisher.graph import BidirectedGraph
from gfapdefinisher.ledger import (
    EndEvidence,
    build_end_evidence_ledger,
    emit_tsv_ledger,
)
from gfapdefinisher.models import SegmentPort, Strand

BENCHMARK_DIR = Path("data/benchmark/DF5120")


def test_end_evidence_dataclass():
    """Verify EndEvidence slots, frozen immutability, and attribute types."""
    port = SegmentPort("utg1", Strand.MINUS)
    ev = EndEvidence(
        segment_name="utg1",
        end="5p",
        port=port,
        length=10000,
        depth=30.0,
        degree=0,
        is_tip=True,
        outgoing_neighbors=(),
        telomere_fwd_count=10,
        telomere_rev_count=0,
        has_terminal_telomere=True,
        rdna_hit_count=0,
        rdna_features=(),
        has_terminal_rdna=False,
        usco_hit_count=1,
        usco_genes=("gene1",),
        tiara_class="eukarya",
    )
    assert ev.segment_name == "utg1"
    assert ev.end == "5p"
    assert ev.port == port
    assert ev.length == 10000
    assert ev.depth == 30.0
    assert ev.degree == 0
    assert ev.is_tip is True
    assert ev.outgoing_neighbors == ()
    assert ev.telomere_fwd_count == 10
    assert ev.telomere_rev_count == 0
    assert ev.has_terminal_telomere is True
    assert ev.rdna_hit_count == 0
    assert ev.rdna_features == ()
    assert ev.has_terminal_rdna is False
    assert ev.usco_hit_count == 1
    assert ev.usco_genes == ("gene1",)
    assert ev.tiara_class == "eukarya"

    with pytest.raises(FrozenInstanceError):
        ev.degree = 1  # type: ignore[misc]


def test_ledger_two_rows_per_unitig_and_deterministic_sorting():
    """Verify that exactly 2 rows (5p and 3p) are emitted per unitig, deterministically sorted."""
    graph = BidirectedGraph()
    graph.add_segment("utgB", length=15000, depth=20.0)
    graph.add_segment("utgA", length=10000, depth=40.0)
    graph.add_segment("utgC", length=25000, depth=10.0)

    ledger = build_end_evidence_ledger(graph, annotations={})

    # Exactly 2 rows per segment
    assert len(ledger) == 6

    # Deterministic sorting by (segment_name, end): utgA 3p, utgA 5p, utgB 3p, utgB 5p, utgC 3p, utgC 5p
    expected_order = [
        ("utgA", "3p"),
        ("utgA", "5p"),
        ("utgB", "3p"),
        ("utgB", "5p"),
        ("utgC", "3p"),
        ("utgC", "5p"),
    ]
    actual_order = [(e.segment_name, e.end) for e in ledger]
    assert actual_order == expected_order

    # Verify ports
    for e in ledger:
        if e.end == "5p":
            assert e.port == SegmentPort(e.segment_name, Strand.MINUS)
        else:
            assert e.port == SegmentPort(e.segment_name, Strand.PLUS)


def test_tip_detection_and_outgoing_neighbors():
    """Verify degree calculation, tip status, and outgoing neighbor names."""
    graph = BidirectedGraph()
    graph.add_segment("utg1", length=5000, depth=30.0)
    graph.add_segment("utg2", length=8000, depth=30.0)
    graph.add_segment("utg3", length=12000, depth=15.0)
    graph.add_segment("utg4", length=12000, depth=15.0)

    # utg1+ -> utg2+ (linear)
    graph.add_edge("utg1+", "utg2+", "50M")
    # utg2+ -> utg3+ and utg2+ -> utg4+ (branching fork)
    graph.add_edge("utg2+", "utg3+", "50M")
    graph.add_edge("utg2+", "utg4+", "50M")

    ledger = build_end_evidence_ledger(graph, annotations={})
    by_key = {(e.segment_name, e.end): e for e in ledger}

    # utg1 5' (utg1-): tip (no incoming into utg1+)
    assert by_key[("utg1", "5p")].degree == 0
    assert by_key[("utg1", "5p")].is_tip is True
    assert by_key[("utg1", "5p")].outgoing_neighbors == ()

    # utg1 3' (utg1+): degree 1 -> utg2
    assert by_key[("utg1", "3p")].degree == 1
    assert by_key[("utg1", "3p")].is_tip is False
    assert by_key[("utg1", "3p")].outgoing_neighbors == ("utg2",)

    # utg2 3' (utg2+): degree 2 -> utg3, utg4
    assert by_key[("utg2", "3p")].degree == 2
    assert by_key[("utg2", "3p")].is_tip is False
    assert set(by_key[("utg2", "3p")].outgoing_neighbors) == {"utg3", "utg4"}

    # utg3 3' (utg3+) and utg4 3' (utg4+): tips
    assert by_key[("utg3", "3p")].degree == 0
    assert by_key[("utg3", "3p")].is_tip is True
    assert by_key[("utg4", "3p")].degree == 0
    assert by_key[("utg4", "3p")].is_tip is True


def test_telomere_and_rdna_end_mapping():
    """Verify telomere windows and rDNA hits map to appropriate 5p vs 3p ends."""
    graph = BidirectedGraph()
    graph.add_segment("utg_arm", length=100000, depth=25.0)

    # Annotations:
    # 5p terminal telomere at 1000, 3p terminal telomere at 95000, internal telomere at 50000
    tidk_wins = [
        TidkWindow("utg_arm", window_start=1000, forward_count=12, reverse_count=2),
        TidkWindow("utg_arm", window_start=50000, forward_count=20, reverse_count=20),
        TidkWindow("utg_arm", window_start=95000, forward_count=0, reverse_count=18),
    ]

    # Barrnap rDNA: 18S near 5p (start=500), 28S near 3p (start=94000), 5S in middle (start=50000)
    barrnap_hits = [
        BarrnapHit("utg_arm", start=500, end=1500, strand="+", feature_name="18S_rRNA", evalue=1e-30),
        BarrnapHit("utg_arm", start=50000, end=50120, strand="+", feature_name="5S_rRNA", evalue=1e-10),
        BarrnapHit("utg_arm", start=94000, end=97000, strand="-", feature_name="28S_rRNA", evalue=1e-40),
    ]

    annotations = {
        "utg_arm": NodeAnnotations(
            sequence_id="utg_arm",
            tidk_windows=tidk_wins,
            barrnap_hits=barrnap_hits,
        )
    }

    # Evaluate with 10 kb terminal window
    ledger = build_end_evidence_ledger(graph, annotations, terminal_window_bp=10000)
    by_end = {e.end: e for e in ledger}

    # 5p end should capture the window at 1000 and the 18S hit at 500
    e_5p = by_end["5p"]
    assert e_5p.telomere_fwd_count == 12
    assert e_5p.telomere_rev_count == 2
    assert e_5p.has_terminal_telomere is True
    assert e_5p.rdna_hit_count == 1
    assert e_5p.rdna_features == ("18S_rRNA",)
    assert e_5p.has_terminal_rdna is True

    # 3p end should capture the window at 95000 and the 28S hit at 94000
    e_3p = by_end["3p"]
    assert e_3p.telomere_fwd_count == 0
    assert e_3p.telomere_rev_count == 18
    assert e_3p.has_terminal_telomere is True
    assert e_3p.rdna_hit_count == 1
    assert e_3p.rdna_features == ("28S_rRNA",)
    assert e_3p.has_terminal_rdna is True


def test_usco_and_tiara_mapping():
    """Verify Compleasm USCO hits and Tiara classifications are attached properly."""
    graph = BidirectedGraph()
    graph.add_segment("utg_gene", length=40000, depth=35.0)

    compleasm_hits = [
        CompleasmHit("USCO_001", "Single", "utg_gene", 5000, 7000, 99.5),
        CompleasmHit("USCO_002", "Single", "utg_gene", 15000, 18000, 98.0),
    ]
    tiara_hit = TiaraHit("utg_gene", "eukarya")

    annotations = {
        "utg_gene": NodeAnnotations(
            sequence_id="utg_gene",
            tiara=tiara_hit,
            compleasm_hits=compleasm_hits,
        )
    }

    ledger = build_end_evidence_ledger(graph, annotations)
    for e in ledger:
        assert e.usco_hit_count == 2
        assert e.usco_genes == ("USCO_001", "USCO_002")
        assert e.tiara_class == "eukarya"


def test_emit_tsv_ledger_formatting(tmp_path: Path):
    """Verify standardized TSV serialization with header and proper columns."""
    graph = BidirectedGraph()
    graph.add_segment("utg1", length=12000, depth=45.5)

    ev_5p = EndEvidence(
        segment_name="utg1",
        end="5p",
        port=SegmentPort("utg1", Strand.MINUS),
        length=12000,
        depth=45.5,
        degree=0,
        is_tip=True,
        outgoing_neighbors=(),
        telomere_fwd_count=15,
        telomere_rev_count=1,
        has_terminal_telomere=True,
        rdna_hit_count=0,
        rdna_features=(),
        has_terminal_rdna=False,
        usco_hit_count=1,
        usco_genes=("geneA",),
        tiara_class="eukarya",
    )
    ev_3p = EndEvidence(
        segment_name="utg1",
        end="3p",
        port=SegmentPort("utg1", Strand.PLUS),
        length=12000,
        depth=45.5,
        degree=2,
        is_tip=False,
        outgoing_neighbors=("utg2", "utg3"),
        telomere_fwd_count=0,
        telomere_rev_count=0,
        has_terminal_telomere=False,
        rdna_hit_count=2,
        rdna_features=("18S_rRNA", "28S_rRNA"),
        has_terminal_rdna=True,
        usco_hit_count=1,
        usco_genes=("geneA",),
        tiara_class="eukarya",
    )

    out_file = tmp_path / "ledger.tsv"
    emit_tsv_ledger([ev_5p, ev_3p], out_file)

    lines = out_file.read_text(encoding="utf-8").splitlines()
    assert len(lines) == 3

    header = lines[0].split("\t")
    expected_cols = [
        "segment_name",
        "end",
        "port",
        "length",
        "depth",
        "degree",
        "is_tip",
        "outgoing_neighbors",
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
    assert header == expected_cols

    row_5p = lines[1].split("\t")
    assert row_5p[0] == "utg1"
    assert row_5p[1] == "5p"
    assert row_5p[2] == "utg1-"
    assert row_5p[3] == "12000"
    assert row_5p[4] == "45.5"
    assert row_5p[5] == "0"
    assert row_5p[6] == "True"
    assert row_5p[7] == ""
    assert row_5p[8] == "15"
    assert row_5p[9] == "1"
    assert row_5p[10] == "True"
    assert row_5p[11] == "0"
    assert row_5p[12] == ""
    assert row_5p[13] == "False"
    assert row_5p[14] == "1"
    assert row_5p[15] == "geneA"
    assert row_5p[16] == "eukarya"

    row_3p = lines[2].split("\t")
    assert row_3p[1] == "3p"
    assert row_3p[2] == "utg1+"
    assert row_3p[5] == "2"
    assert row_3p[6] == "False"
    assert row_3p[7] == "utg2,utg3"
    assert row_3p[12] == "18S_rRNA,28S_rRNA"

    # Also test stream writing with io.StringIO
    buf = io.StringIO()
    emit_tsv_ledger([ev_5p], buf)
    buf_lines = buf.getvalue().splitlines()
    assert len(buf_lines) == 2


def test_cli_report_subcommand(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    """Verify gfapdefinisher report CLI runs and emits output file or prints to stdout."""
    # Write synthetic GFA
    gfa_path = tmp_path / "test.gfa"
    gfa_content = (
        "H\tVN:Z:1.0\n"
        "S\tutg001\tACGTACGT\tLN:i:8\trd:f:20.5\n"
        "S\tutg002\tTGCATGCA\tLN:i:8\trd:f:10.0\n"
        "L\tutg001\t+\tutg002\t+\t0M\n"
    )
    gfa_path.write_text(gfa_content, encoding="utf-8")

    # Write synthetic Tiara
    tiara_path = tmp_path / "tiara.txt"
    tiara_path.write_text("sequence_id\tclass_fst_stage\nutg001\teukarya\n", encoding="utf-8")

    # Write synthetic TIDK
    tidk_path = tmp_path / "tidk.tsv"
    tidk_path.write_text(
        "id\twindow\tforward_repeat_number\treverse_repeat_number\nutg001\t0\t5\t0\n", encoding="utf-8"
    )

    out_tsv = tmp_path / "output_ledger.tsv"

    # Run with file output
    cmd = [
        "report",
        "--gfa",
        str(gfa_path),
        "--tiara",
        str(tiara_path),
        "--tidk",
        str(tidk_path),
        "--out",
        str(out_tsv),
    ]
    main(cmd)

    assert out_tsv.exists()
    content = out_tsv.read_text(encoding="utf-8")
    assert "utg001" in content
    assert "utg002" in content
    assert "eukarya" in content

    # Run with stdout output
    main(["report", "--gfa", str(gfa_path)])
    captured = capsys.readouterr()
    assert "segment_name\tend\tport" in captured.out
    assert "utg001" in captured.out


@pytest.mark.benchmark
def test_benchmark_df5120_staged_integration(tmp_path: Path):
    """Integration test running against staged benchmark files if present."""
    if not BENCHMARK_DIR.exists():
        pytest.skip(f"Benchmark directory {BENCHMARK_DIR} not present.")

    gfa_candidates = [
        BENCHMARK_DIR / "asm.p_utg.gfa.gz",
        BENCHMARK_DIR / "asm.p_utg.gfa",
    ]
    gfa_path = next((p for p in gfa_candidates if p.exists()), None)
    if gfa_path is None:
        pytest.skip("No GFA file found in benchmark directory.")

    tiara_path = BENCHMARK_DIR / "tiara_unitigs.txt"
    tidk_path = BENCHMARK_DIR / "tidk_telomeres_telomeric_repeat_windows.tsv"
    barrnap_path = BENCHMARK_DIR / "barrnap_rdna.gff"
    compleasm_path = BENCHMARK_DIR / "compleasm_unitigs_full_table.tsv"

    out_tsv = tmp_path / "benchmark_ledger.tsv"
    cmd = [
        "report",
        "--gfa",
        str(gfa_path),
        "--out",
        str(out_tsv),
    ]
    if tiara_path.exists():
        cmd.extend(["--tiara", str(tiara_path)])
    if tidk_path.exists():
        cmd.extend(["--tidk", str(tidk_path)])
    if barrnap_path.exists():
        cmd.extend(["--barrnap", str(barrnap_path)])
    if compleasm_path.exists():
        cmd.extend(["--compleasm", str(compleasm_path)])

    main(cmd)

    assert out_tsv.exists()
    lines = out_tsv.read_text(encoding="utf-8").splitlines()
    assert len(lines) > 10
    header = lines[0].split("\t")
    assert "segment_name" in header
    assert "has_terminal_telomere" in header
