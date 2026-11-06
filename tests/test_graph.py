"""Tests for BidirectedGraph and lossless GFA1 loading."""

import gzip
from pathlib import Path

import pytest

from gfapdefinisher.graph import BidirectedGraph, load_gfa
from gfapdefinisher.models import SegmentPort, Strand


def test_graph_initialization_and_segment_addition():
    """Verify segment addition, attributes, and indexing."""
    graph = BidirectedGraph()
    seg = graph.add_segment("utg001", length=15000, depth=32.5, sequence="A" * 15000)

    assert seg.name == "utg001"
    assert seg.length == 15000
    assert seg.depth == 32.5
    assert seg.sequence == "A" * 15000
    assert len(graph) == 1
    assert "utg001" in graph
    assert graph["utg001"] == seg

    # Segments with tags
    graph.add_segment("utg002", length=20000, depth=15.0, tags={"LN": 20000, "rd": 15, "HG": "p"})
    assert graph.tags["utg002"]["HG"] == "p"


def test_linear_overlap_topology_and_duality():
    """Verify linear overlap topology (A -> B) and Watson-Crick duality."""
    graph = BidirectedGraph()
    graph.add_segment("A", length=1000)
    graph.add_segment("B", length=2000)

    p_a_fwd = SegmentPort("A", Strand.FORWARD)
    p_b_fwd = SegmentPort("B", Strand.FORWARD)

    edge = graph.add_edge(p_a_fwd, p_b_fwd, cigar="100M", l1=1000, l2=2000)

    assert edge.from_port == p_a_fwd
    assert edge.to_port == p_b_fwd
    assert edge.cigar == "100M"
    assert edge.overlap_len == 100

    # Watson-Crick duality check:
    # get_outgoing(u, +) matches get_incoming(v, +)
    outgoing_a_fwd = graph.get_outgoing("A", "+")
    incoming_b_fwd = graph.get_incoming("B", "+")
    assert outgoing_a_fwd == [edge]
    assert incoming_b_fwd == [edge]

    # Reverse dual check:
    # get_outgoing(v, -) matches get_incoming(u, -)
    outgoing_b_rev = graph.get_outgoing("B", "-")
    incoming_a_rev = graph.get_incoming("A", "-")
    assert len(outgoing_b_rev) == 1
    assert len(incoming_a_rev) == 1
    assert outgoing_b_rev[0] == edge.dual
    assert incoming_a_rev[0] == edge.dual
    assert outgoing_b_rev == incoming_a_rev

    # Tip queries
    # A- is a 5' contig tip; B+ is a 3' contig tip
    assert graph.is_tip("A", "-") is True
    assert graph.is_tip("B", "+") is True
    assert graph.is_tip("A", "+") is False
    assert graph.is_tip("B", "-") is False

    # Degrees
    assert graph.degree("A", "+") == 1
    assert graph.degree("A", "-") == 0
    assert graph.degree("B", "+") == 0
    assert graph.degree("B", "-") == 1


def test_bubble_topology():
    """Verify bubble topology: A -> B -> D and A -> C -> D."""
    graph = BidirectedGraph()
    for name in ["A", "B", "C", "D"]:
        graph.add_segment(name, length=500)

    # A forks into B and C
    graph.add_edge("A+", "B+", "50M")
    graph.add_edge("A+", "C+", "50M")

    # B and C merge into D
    graph.add_edge("B+", "D+", "50M")
    graph.add_edge("C+", "D+", "50M")

    # Outgoing degree of A+ is 2 (branching)
    assert graph.degree("A", "+") == 2
    assert graph.is_tip("A", "+") is False
    out_a = {e.to_port.name for e in graph.get_outgoing("A", "+")}
    assert out_a == {"B", "C"}

    # B+ and C+ have degree 1 (connecting to D+)
    assert graph.degree("B", "+") == 1
    assert graph.degree("C", "+") == 1

    # D+ is a tip (sink)
    assert graph.degree("D", "+") == 0
    assert graph.is_tip("D", "+") is True

    # Incoming to D+ has 2 edges (from B+ and C+)
    inc_d = {e.from_port.name for e in graph.get_incoming("D", "+")}
    assert inc_d == {"B", "C"}

    # Dual bubble in reverse orientation: D- forks into B- and C-, which merge into A-
    assert graph.degree("D", "-") == 2
    out_d_rev = {e.to_port.name for e in graph.get_outgoing("D", "-")}
    assert out_d_rev == {"B", "C"}
    assert graph.degree("A", "-") == 0
    assert graph.is_tip("A", "-") is True


def test_cycle_topology():
    """Verify cycle topology: A -> B -> C -> A."""
    graph = BidirectedGraph()
    for name in ["A", "B", "C"]:
        graph.add_segment(name, length=300)

    graph.add_edge("A+", "B+", "30M")
    graph.add_edge("B+", "C+", "30M")
    graph.add_edge("C+", "A+", "30M")

    # No tips should exist in a closed cycle
    for name in ["A", "B", "C"]:
        assert graph.degree(name, "+") == 1
        assert graph.degree(name, "-") == 1
        assert graph.is_tip(name, "+") is False
        assert graph.is_tip(name, "-") is False


def test_fork_topology():
    """Verify fork topology: A -> B and A -> C (divergence without reconvergence)."""
    graph = BidirectedGraph()
    for name in ["A", "B", "C"]:
        graph.add_segment(name, length=400)

    graph.add_edge("A+", "B+", "40M")
    graph.add_edge("A+", "C+", "40M")

    assert graph.degree("A", "+") == 2
    assert graph.degree("B", "+") == 0
    assert graph.degree("C", "+") == 0

    assert graph.is_tip("B", "+") is True
    assert graph.is_tip("C", "+") is True
    assert graph.is_tip("A", "+") is False


def test_load_gfa_with_and_without_sequences(tmp_path: Path):
    """Verify GFA loading with sequences loaded or holding metadata only."""
    gfa_content = "\n".join(
        [
            "H\tVN:Z:1.0",
            "S\tutg001\tACGTACGT\tLN:i:8\trd:i:25\tHG:A:m",
            "S\tutg002\tTTTTAAAA\tLN:i:8\trd:i:30\tHG:A:p",
            "L\tutg001\t+\tutg002\t+\t4M\tL1:i:8\tL2:i:8",
            "A\tutg001\t0\t+\tread_01\t0\t1000\t*",
        ]
    )
    gfa_file = tmp_path / "test.gfa"
    gfa_file.write_text(gfa_content, encoding="utf-8")

    # 1. Load metadata only (load_sequences=False)
    graph_meta = load_gfa(gfa_file, load_sequences=False)
    assert len(graph_meta) == 2
    assert graph_meta["utg001"].sequence is None
    assert graph_meta["utg001"].length == 8
    assert graph_meta["utg001"].depth == 25.0
    assert graph_meta.tags["utg001"]["HG"] == "m"
    assert len(graph_meta.a_lines) == 1
    assert "read_01" in graph_meta.a_lines[0]

    # Lazy on-demand loading
    seq1 = graph_meta.get_sequence("utg001")
    assert seq1 == "ACGTACGT"
    assert graph_meta["utg001"].sequence == "ACGTACGT"

    # 2. Load with sequences directly
    graph_full = load_gfa(gfa_file, load_sequences=True)
    assert graph_full["utg001"].sequence == "ACGTACGT"
    assert graph_full["utg002"].sequence == "TTTTAAAA"


def test_load_gfa_gzip_support(tmp_path: Path):
    """Verify loading from compressed GFA (.gfa.gz)."""
    gfa_content = "\n".join(
        [
            "H\tVN:Z:1.0",
            "S\tseg1\t*\tLN:i:1000\tdp:f:18.5",
            "S\tseg2\t*\tLN:i:2000\tdp:f:21.0",
            "L\tseg1\t+\tseg2\t+\t50M",
        ]
    )
    gz_file = tmp_path / "compressed.gfa.gz"
    with gzip.open(gz_file, "wt", encoding="utf-8") as f:
        f.write(gfa_content)

    graph = load_gfa(gz_file)
    assert len(graph) == 2
    assert graph["seg1"].length == 1000
    assert graph["seg1"].depth == 18.5
    assert graph["seg2"].length == 2000
    assert graph.degree("seg1", "+") == 1
    assert graph.degree("seg2", "+") == 0


def test_gfa_roundtrip_lossless_preservation(tmp_path: Path):
    """Verify serialization and roundtrip preservation of all tags and A lines."""
    gfa_content = "\n".join(
        [
            "H\tVN:Z:1.0",
            "S\tutg0001\tACGT\tLN:i:4\trd:i:40\tHG:A:p\tFC:i:1",
            "S\tutg0002\tTGCA\tLN:i:4\trd:i:35\tHG:A:m\tFC:i:2",
            "L\tutg0001\t+\tutg0002\t+\t2M\tL1:i:4\tL2:i:4",
            "A\tutg0001\t0\t+\tr1\t0\t100\t*",
        ]
    )
    src_file = tmp_path / "original.gfa"
    src_file.write_text(gfa_content, encoding="utf-8")

    graph = load_gfa(src_file, load_sequences=True)
    out_file = tmp_path / "roundtrip.gfa"
    graph.write_gfa(out_file)

    reloaded = load_gfa(out_file, load_sequences=True)
    assert len(reloaded) == len(graph)
    assert reloaded["utg0001"].sequence == graph["utg0001"].sequence
    assert reloaded["utg0001"].length == graph["utg0001"].length
    assert reloaded["utg0001"].depth == graph["utg0001"].depth
    assert reloaded.tags["utg0001"]["HG"] == "p"
    assert reloaded.tags["utg0001"]["FC"] == 1
    assert reloaded.tags["utg0002"]["FC"] == 2
    assert len(reloaded.a_lines) == 1
    assert len(reloaded.canonical_edges) == 1


def test_cigar_with_indels_in_graph_edge():
    """Verify dovetail edge computation with complex CIGAR indels."""
    graph = BidirectedGraph()
    graph.add_segment("u1", length=500)
    graph.add_segment("u2", length=500)

    # Overlap with insertion
    edge = graph.add_edge("u1+", "u2+", "10M2I10M")
    # cigar_overlap_length for 10M2I10M: target consumption=20, query consumption=22, max=22
    assert edge.overlap_len == 22
    assert graph.degree("u1", "+") == 1
    assert graph.degree("u2", "-") == 1


def test_load_gfa_reciprocal_dual_deduplication(tmp_path: Path):
    """Verify that reciprocal dual L lines in GFA do not inflate degree or duplicate edges."""
    gfa_content = "\n".join(
        [
            "H\tVN:Z:1.0",
            "S\tsegA\t*\tLN:i:1000",
            "S\tsegB\t*\tLN:i:2000",
            "S\tsegC\t*\tLN:i:1500",
            "L\tsegA\t+\tsegB\t+\t50M",
            "L\tsegB\t-\tsegA\t-\t50M",
            "L\tsegB\t+\tsegC\t+\t40M",
            "L\tsegC\t-\tsegB\t-\t40M",
        ]
    )
    gfa_file = tmp_path / "reciprocal.gfa"
    gfa_file.write_text(gfa_content, encoding="utf-8")

    graph = load_gfa(gfa_file)

    # Exactly 2 canonical edges
    assert len(graph.canonical_edges) == 2

    # Port outgoing degrees must NOT be doubled (no inflation):
    assert graph.degree("segA", "+") == 1
    assert graph.degree("segA", "-") == 0  # 5' tip
    assert graph.degree("segB", "+") == 1
    assert graph.degree("segB", "-") == 1
    assert graph.degree("segC", "+") == 0  # 3' tip
    assert graph.degree("segC", "-") == 1

    # Incoming degrees must not be doubled:
    assert len(graph.get_incoming("segB", "+")) == 1
    assert len(graph.get_incoming("segC", "+")) == 1
    assert len(graph.get_incoming("segA", "-")) == 1
    assert len(graph.get_incoming("segB", "-")) == 1

    # Outgoing edge count per port:
    assert len(graph.get_outgoing("segA", "+")) == 1
    assert len(graph.get_outgoing("segB", "-")) == 1


def test_forward_only_vs_reciprocal_gfa_equivalence(tmp_path: Path):
    """Verify graphs loaded from forward-only vs reciprocal GFA have identical topology and degrees."""
    fwd_gfa = "\n".join(
        [
            "H\tVN:Z:1.0",
            "S\tu1\t*\tLN:i:500",
            "S\tu2\t*\tLN:i:600",
            "L\tu1\t+\tu2\t+\t30M",
        ]
    )
    recip_gfa = "\n".join(
        [
            "H\tVN:Z:1.0",
            "S\tu1\t*\tLN:i:500",
            "S\tu2\t*\tLN:i:600",
            "L\tu1\t+\tu2\t+\t30M",
            "L\tu2\t-\tu1\t-\t30M",
        ]
    )
    p_fwd = tmp_path / "fwd.gfa"
    p_recip = tmp_path / "recip.gfa"
    p_fwd.write_text(fwd_gfa, encoding="utf-8")
    p_recip.write_text(recip_gfa, encoding="utf-8")

    g_fwd = load_gfa(p_fwd)
    g_recip = load_gfa(p_recip)

    for seg in ["u1", "u2"]:
        for strand in ["+", "-"]:
            assert g_fwd.degree(seg, strand) == g_recip.degree(seg, strand)
            assert g_fwd.is_tip(seg, strand) == g_recip.is_tip(seg, strand)

    assert len(g_fwd.canonical_edges) == len(g_recip.canonical_edges) == 1


def test_coordinate_l2_zero_preserved():
    """Verify that l2 = 0 is preserved and not overwritten with segment length."""
    graph = BidirectedGraph()
    graph.add_segment("s1", length=1000)
    graph.add_segment("s2", length=2000)

    edge = graph.add_edge("s1+", "s2+", cigar="50M", l1=1000, l2=0)
    assert edge.l1 == 1000
    assert edge.l2 == 0


def test_edge_cigar_upgrade_and_metadata_merge():
    """Verify adding edge with better CIGAR upgrades placeholder edge and merges tags."""
    graph = BidirectedGraph()
    graph.add_segment("s1", length=1000)
    graph.add_segment("s2", length=2000)

    # First add with placeholder CIGAR
    edge1 = graph.add_edge("s1+", "s2+", cigar="*", tags={"t1": 1})
    assert edge1.cigar == "*"

    # Now add with specific CIGAR and additional tag
    edge2 = graph.add_edge("s1+", "s2+", cigar="50M", tags={"t2": 2})
    assert edge2.cigar == "50M"
    assert edge2.overlap_len == 50

    # Outgoing degree remains 1
    assert graph.degree("s1", "+") == 1
    out_edge = graph.get_outgoing("s1", "+")[0]
    assert out_edge.cigar == "50M"
    assert out_edge.overlap_len == 50

    # Tags merged
    tags = graph.edge_tags[(SegmentPort("s1", Strand.FORWARD), SegmentPort("s2", Strand.FORWARD))]
    assert tags["t1"] == 1
    assert tags["t2"] == 2


@pytest.mark.benchmark
def test_benchmark_df5120_degree_distribution():
    """Verify exact degree distribution on benchmark DF5120 assembly graph."""
    bench_candidates = [
        Path("data/benchmark/DF5120/asm.p_utg.gfa.gz"),
        Path("data/benchmark/DF5120/asm.p_utg.gfa"),
    ]
    bench_file = next((p for p in bench_candidates if p.exists()), None)
    if bench_file is None:
        pytest.skip("Benchmark DF5120 GFA file not present.")

    graph = load_gfa(bench_file)

    # Port outgoing degree distribution in raw GFA:
    # deg 0 (tips): 233 ports
    # deg 1 (collinear): 141 ports
    # deg 2 (forks): 69 ports
    # deg 3 (complex branches): 3 ports
    degree_counts: dict[int, int] = {}
    for seg_name in graph.segments:
        for strand in [Strand.FORWARD, Strand.REVERSE]:
            deg = graph.degree(seg_name, strand)
            degree_counts[deg] = degree_counts.get(deg, 0) + 1

    assert degree_counts.get(0, 0) == 233
    assert degree_counts.get(1, 0) == 141
    assert degree_counts.get(2, 0) == 69
    assert degree_counts.get(3, 0) == 3
    assert len(graph.canonical_edges) == 144
