"""Tests for rDNA cassette graph-walking resolver and head-to-tail copy chaining."""

from io import StringIO
from pathlib import Path

import pytest

from gfapdefinisher.annotations import BarrnapHit, NodeAnnotations, load_barrnap
from gfapdefinisher.graph import BidirectedGraph, load_gfa
from gfapdefinisher.ledger import EndEvidence
from gfapdefinisher.models import SegmentPort, Strand
from gfapdefinisher.rdna import RdnaCassetteResolver, is_rdna_unitig
from gfapdefinisher.resolutions import (
    EndState,
    IntervalQuirk,
    parse_end_states_tsv,
    parse_quirks_tsv,
)

BENCHMARK_DIR = Path("data/benchmark/DF5120")


def _make_evidence(
    segment_name: str,
    end: str,
    depth: float = 60.0,
    length: int = 10000,
    has_terminal_telomere: bool = False,
    rdna_hit_count: int = 0,
    has_terminal_rdna: bool = False,
    rdna_features: tuple[str, ...] = (),
    usco_hit_count: int = 0,
) -> EndEvidence:
    strand = Strand.MINUS if end == "5p" else Strand.PLUS
    port = SegmentPort(segment_name, strand)
    return EndEvidence(
        segment_name=segment_name,
        end=end,
        port=port,
        length=length,
        depth=depth,
        degree=0,
        is_tip=True,
        outgoing_neighbors=(),
        has_terminal_telomere=has_terminal_telomere,
        rdna_hit_count=rdna_hit_count,
        rdna_features=rdna_features,
        has_terminal_rdna=has_terminal_rdna,
        usco_hit_count=usco_hit_count,
        tiara_class="eukarya",
    )


def test_find_rdna_entry_ports_identifies_core_stems():
    """Verify core stems connecting into rDNA clusters are correctly identified as entry ports."""
    graph = BidirectedGraph()
    graph.add_segment("chr_stem", length=50000, depth=60.0)
    graph.add_segment("rdna_unitig", length=15000, depth=240.0)
    graph.add_segment("bubble_unitig", length=10000, depth=60.0)

    # Stem connects to rDNA and bubble
    p_stem = SegmentPort("chr_stem", Strand.FORWARD)
    p_rdna = SegmentPort("rdna_unitig", Strand.FORWARD)
    p_bubble = SegmentPort("bubble_unitig", Strand.FORWARD)

    graph.add_edge(p_stem, p_rdna, cigar="100M")
    graph.add_edge(p_stem, p_bubble, cigar="100M")

    ev_stem_5p = _make_evidence("chr_stem", "5p")
    ev_stem_3p = _make_evidence("chr_stem", "3p")
    ev_rdna_5p = _make_evidence("rdna_unitig", "5p", rdna_hit_count=3, has_terminal_rdna=True)
    ev_rdna_3p = _make_evidence("rdna_unitig", "3p", rdna_hit_count=3, has_terminal_rdna=True)
    ev_bub_5p = _make_evidence("bubble_unitig", "5p")
    ev_bub_3p = _make_evidence("bubble_unitig", "3p")

    ledger = {
        "chr_stem+": ev_stem_3p,
        "chr_stem-": ev_stem_5p,
        "rdna_unitig+": ev_rdna_3p,
        "rdna_unitig-": ev_rdna_5p,
        "bubble_unitig+": ev_bub_3p,
        "bubble_unitig-": ev_bub_5p,
    }

    resolver = RdnaCassetteResolver(graph, ledger)
    entry_ports = resolver.find_rdna_entry_ports()

    assert len(entry_ports) == 1
    assert entry_ports[0] == p_stem
    assert resolver.identify_entry_anchor() == p_stem


def test_head_to_tail_traversal_and_copy_chaining():
    """Verify traversal chains multiple rDNA unitigs with consistent 18S->5.8S->28S orientation."""
    graph = BidirectedGraph()
    graph.add_segment("stem", length=20000, depth=60.0, sequence="A" * 20000)
    graph.add_segment("rdna_1", length=10000, depth=200.0, sequence="C" * 10000)
    graph.add_segment("rdna_2", length=10000, depth=200.0, sequence="G" * 10000)
    graph.add_segment("rdna_3", length=10000, depth=200.0, sequence="C" * 10000)
    graph.add_segment("distal_cap", length=5000, depth=30.0, sequence="T" * 5000)

    p_stem = SegmentPort("stem", Strand.FORWARD)
    p_r1 = SegmentPort("rdna_1", Strand.FORWARD)
    p_r2 = SegmentPort("rdna_2", Strand.FORWARD)
    p_r3 = SegmentPort("rdna_3", Strand.FORWARD)
    p_cap = SegmentPort("distal_cap", Strand.FORWARD)

    graph.add_edge(p_stem, p_r1, cigar="100M")
    graph.add_edge(p_r1, p_r2, cigar="100M")
    graph.add_edge(p_r2, p_r3, cigar="100M")
    graph.add_edge(p_r3, p_cap, cigar="100M")

    # Annotations establishing intrinsic forward orientation for all rDNA copies
    annotations = {
        "rdna_1": NodeAnnotations(
            sequence_id="rdna_1",
            barrnap_hits=[
                BarrnapHit("rdna_1", 100, 1800, "+", "18S_rRNA", 1e-40),
                BarrnapHit("rdna_1", 2000, 2150, "+", "5_8S_rRNA", 1e-10),
                BarrnapHit("rdna_1", 2500, 6000, "+", "28S_rRNA", 1e-60),
            ],
        ),
        "rdna_2": NodeAnnotations(
            sequence_id="rdna_2",
            barrnap_hits=[
                BarrnapHit("rdna_2", 100, 1800, "+", "18S_rRNA", 1e-40),
                BarrnapHit("rdna_2", 2500, 6000, "+", "28S_rRNA", 1e-60),
            ],
        ),
        "rdna_3": NodeAnnotations(
            sequence_id="rdna_3",
            barrnap_hits=[
                BarrnapHit("rdna_3", 100, 1800, "+", "18S_rRNA", 1e-40),
                BarrnapHit("rdna_3", 2500, 6000, "+", "28S_rRNA", 1e-60),
            ],
        ),
    }

    ledger = {
        "stem+": _make_evidence("stem", "3p"),
        "distal_cap+": _make_evidence("distal_cap", "3p", has_terminal_telomere=True),
    }

    resolver = RdnaCassetteResolver(graph, ledger, annotations=annotations)

    # 1. Test traverse_rdna_cassette chains rDNA unitigs
    rdna_chain = resolver.traverse_rdna_cassette(p_stem)
    assert rdna_chain == [p_r1, p_r2, p_r3]

    # 2. Test reach_distal_telomere_cap identifies full path
    full_path = resolver.reach_distal_telomere_cap(p_stem)
    assert full_path == [p_stem, p_r1, p_r2, p_r3, p_cap]

    # 3. Test assemble_contiguous_path
    assembled = resolver.assemble_contiguous_path(full_path)
    assert assembled is not None
    assert len(assembled) == 20000 + (10000 - 100) * 3 + (5000 - 100)


def test_inversion_switch_rejection():
    """Verify that inverted rDNA unitigs (28S->18S) are rejected during graph walking."""
    graph = BidirectedGraph()
    graph.add_segment("stem", length=20000, depth=60.0)
    graph.add_segment("rdna_fwd", length=10000, depth=200.0)
    graph.add_segment("rdna_inv", length=10000, depth=200.0)
    graph.add_segment("cap_inv", length=5000, depth=30.0)
    graph.add_segment("rdna_valid", length=10000, depth=200.0)
    graph.add_segment("cap_valid", length=5000, depth=30.0)

    p_stem = SegmentPort("stem", Strand.FORWARD)
    p_fwd = SegmentPort("rdna_fwd", Strand.FORWARD)
    p_inv = SegmentPort("rdna_inv", Strand.REVERSE)  # Inverted port
    p_cap_inv = SegmentPort("cap_inv", Strand.FORWARD)
    p_valid = SegmentPort("rdna_valid", Strand.FORWARD)
    p_cap_valid = SegmentPort("cap_valid", Strand.FORWARD)

    # Stem -> rdna_fwd
    graph.add_edge(p_stem, p_fwd, cigar="100M")

    # Branch A: rdna_fwd -> rdna_inv- (inversion switch: rdna_inv intrinsic is +)
    graph.add_edge(p_fwd, p_inv, cigar="100M")
    graph.add_edge(p_inv, p_cap_inv, cigar="100M")

    # Branch B: rdna_fwd -> rdna_valid+ -> cap_valid+
    graph.add_edge(p_fwd, p_valid, cigar="100M")
    graph.add_edge(p_valid, p_cap_valid, cigar="100M")

    # rdna_inv has intrinsic forward orientation, so traversing as reverse is an inversion switch!
    annotations = {
        "rdna_fwd": NodeAnnotations(
            sequence_id="rdna_fwd",
            barrnap_hits=[
                BarrnapHit("rdna_fwd", 100, 1800, "+", "18S_rRNA", 1e-40),
                BarrnapHit("rdna_fwd", 2500, 6000, "+", "28S_rRNA", 1e-60),
            ],
        ),
        "rdna_inv": NodeAnnotations(
            sequence_id="rdna_inv",
            barrnap_hits=[
                BarrnapHit("rdna_inv", 100, 1800, "+", "18S_rRNA", 1e-40),
                BarrnapHit("rdna_inv", 2500, 6000, "+", "28S_rRNA", 1e-60),
            ],
        ),
        "rdna_valid": NodeAnnotations(
            sequence_id="rdna_valid",
            barrnap_hits=[
                BarrnapHit("rdna_valid", 100, 1800, "+", "18S_rRNA", 1e-40),
                BarrnapHit("rdna_valid", 2500, 6000, "+", "28S_rRNA", 1e-60),
            ],
        ),
    }

    ledger = {
        "stem+": _make_evidence("stem", "3p"),
        "cap_inv+": _make_evidence("cap_inv", "3p", has_terminal_telomere=True),
        "cap_valid+": _make_evidence("cap_valid", "3p", has_terminal_telomere=True),
    }

    resolver = RdnaCassetteResolver(graph, ledger, annotations=annotations)

    # Inversion check
    assert resolver.is_consistent_rdna_orientation(p_inv) is False
    assert resolver.is_consistent_rdna_orientation(p_valid) is True

    # Traversal must reject Branch A and choose Branch B
    full_path = resolver.reach_distal_telomere_cap(p_stem)
    assert full_path is not None
    assert full_path == [p_stem, p_fwd, p_valid, p_cap_valid]


def test_cycle_prevention_no_infinite_loop():
    """Verify that cyclic repeat loops in rDNA tangles are terminated by visited tracking."""
    graph = BidirectedGraph()
    graph.add_segment("stem", length=20000, depth=60.0)
    graph.add_segment("rdna_cycle_1", length=10000, depth=300.0)
    graph.add_segment("rdna_cycle_2", length=10000, depth=300.0)
    graph.add_segment("cap", length=5000, depth=30.0)

    p_stem = SegmentPort("stem", Strand.FORWARD)
    p_c1 = SegmentPort("rdna_cycle_1", Strand.FORWARD)
    p_c2 = SegmentPort("rdna_cycle_2", Strand.FORWARD)
    p_cap = SegmentPort("cap", Strand.FORWARD)

    # Cycle: c1 -> c2 -> c1
    graph.add_edge(p_stem, p_c1, cigar="100M")
    graph.add_edge(p_c1, p_c2, cigar="100M")
    graph.add_edge(p_c2, p_c1, cigar="100M")
    # Exit to cap
    graph.add_edge(p_c2, p_cap, cigar="100M")

    annotations = {
        "rdna_cycle_1": NodeAnnotations(
            sequence_id="rdna_cycle_1",
            barrnap_hits=[BarrnapHit("rdna_cycle_1", 100, 1800, "+", "18S_rRNA", 1e-40)],
        ),
        "rdna_cycle_2": NodeAnnotations(
            sequence_id="rdna_cycle_2",
            barrnap_hits=[BarrnapHit("rdna_cycle_2", 100, 1800, "+", "18S_rRNA", 1e-40)],
        ),
    }

    ledger = {
        "stem+": _make_evidence("stem", "3p"),
        "cap+": _make_evidence("cap", "3p", has_terminal_telomere=True),
    }

    resolver = RdnaCassetteResolver(graph, ledger, annotations=annotations)
    full_path = resolver.reach_distal_telomere_cap(p_stem)

    # Path must reach cap without duplicating segments
    assert full_path is not None
    assert full_path == [p_stem, p_c1, p_c2, p_cap]


def test_distal_cap_reachability_and_selection_prefers_longest_chain():
    """Verify that resolver prefers path chaining as many consistent rDNA unitigs as possible."""
    graph = BidirectedGraph()
    graph.add_segment("stem", length=20000, depth=60.0)
    graph.add_segment("rdna_short", length=10000, depth=200.0)
    graph.add_segment("rdna_1", length=10000, depth=200.0)
    graph.add_segment("rdna_2", length=10000, depth=200.0)
    graph.add_segment("rdna_3", length=10000, depth=200.0)
    graph.add_segment("cap1", length=5000, depth=30.0)
    graph.add_segment("cap2", length=5000, depth=30.0)

    p_stem = SegmentPort("stem", Strand.FORWARD)
    p_short = SegmentPort("rdna_short", Strand.FORWARD)
    p_r1 = SegmentPort("rdna_1", Strand.FORWARD)
    p_r2 = SegmentPort("rdna_2", Strand.FORWARD)
    p_r3 = SegmentPort("rdna_3", Strand.FORWARD)
    p_cap1 = SegmentPort("cap1", Strand.FORWARD)
    p_cap2 = SegmentPort("cap2", Strand.FORWARD)

    # Short path: stem -> rdna_short -> cap1 (1 copy)
    graph.add_edge(p_stem, p_short, cigar="100M")
    graph.add_edge(p_short, p_cap1, cigar="100M")

    # Long path: stem -> rdna_1 -> rdna_2 -> rdna_3 -> cap2 (3 copies)
    graph.add_edge(p_stem, p_r1, cigar="100M")
    graph.add_edge(p_r1, p_r2, cigar="100M")
    graph.add_edge(p_r2, p_r3, cigar="100M")
    graph.add_edge(p_r3, p_cap2, cigar="100M")

    annotations = {
        "rdna_short": NodeAnnotations(
            sequence_id="rdna_short",
            barrnap_hits=[BarrnapHit("rdna_short", 100, 1800, "+", "18S_rRNA", 1e-40)],
        ),
        "rdna_1": NodeAnnotations(
            sequence_id="rdna_1",
            barrnap_hits=[BarrnapHit("rdna_1", 100, 1800, "+", "18S_rRNA", 1e-40)],
        ),
        "rdna_2": NodeAnnotations(
            sequence_id="rdna_2",
            barrnap_hits=[BarrnapHit("rdna_2", 100, 1800, "+", "18S_rRNA", 1e-40)],
        ),
        "rdna_3": NodeAnnotations(
            sequence_id="rdna_3",
            barrnap_hits=[BarrnapHit("rdna_3", 100, 1800, "+", "18S_rRNA", 1e-40)],
        ),
    }

    ledger = {
        "stem+": _make_evidence("stem", "3p"),
        "cap1+": _make_evidence("cap1", "3p", has_terminal_telomere=True),
        "cap2+": _make_evidence("cap2", "3p", has_terminal_telomere=True),
    }

    resolver = RdnaCassetteResolver(graph, ledger, annotations=annotations)
    full_path = resolver.reach_distal_telomere_cap(p_stem)

    # Must prefer 3-copy path over 1-copy path
    assert full_path is not None
    assert full_path == [p_stem, p_r1, p_r2, p_r3, p_cap2]


def test_unreached_distal_cap_resolution_is_unresolved():
    """Verify BIO-001: naked rDNA ends without telomeric cap remain unresolved."""
    graph = BidirectedGraph()
    graph.add_segment("stem", length=20000, depth=60.0)
    graph.add_segment("rdna_end", length=10000, depth=200.0)

    p_stem = SegmentPort("stem", Strand.FORWARD)
    p_rdna = SegmentPort("rdna_end", Strand.FORWARD)
    graph.add_edge(p_stem, p_rdna, cigar="100M")

    annotations = {
        "rdna_end": NodeAnnotations(
            sequence_id="rdna_end",
            barrnap_hits=[BarrnapHit("rdna_end", 100, 1800, "+", "18S_rRNA", 1e-40)],
        ),
    }

    ledger = {
        "stem+": _make_evidence("stem", "3p"),
        "rdna_end+": _make_evidence("rdna_end", "3p", rdna_hit_count=1),
    }

    resolver = RdnaCassetteResolver(graph, ledger, annotations=annotations)
    assert resolver.reach_distal_telomere_cap(p_stem) is None

    res = resolver.resolve_rdna_cassette(chrom="II", side="5p", entry_port=p_stem)
    assert res.state == EndState.END_SHARED_EXEMPT
    assert res.end_state_record.resolved is False
    assert res.quirk_record.quirk_type == IntervalQuirk.RDNA_CASSETTE
    assert res.quirk_record.right_anchor_port is None


def test_tsv_emission_and_serde_roundtrip():
    """Verify TSV serialization and deserialization of EndStateRecord and QuirkRecord."""
    graph = BidirectedGraph()
    graph.add_segment("stem", length=20000, depth=60.0)
    graph.add_segment("rdna_1", length=10000, depth=200.0)
    graph.add_segment("cap", length=5000, depth=30.0)

    p_stem = SegmentPort("stem", Strand.FORWARD)
    p_r1 = SegmentPort("rdna_1", Strand.FORWARD)
    p_cap = SegmentPort("cap", Strand.FORWARD)

    graph.add_edge(p_stem, p_r1, cigar="100M")
    graph.add_edge(p_r1, p_cap, cigar="100M")

    annotations = {
        "rdna_1": NodeAnnotations(
            sequence_id="rdna_1",
            barrnap_hits=[BarrnapHit("rdna_1", 100, 1800, "+", "18S_rRNA", 1e-40)],
        ),
    }

    ledger = {
        "stem+": _make_evidence("stem", "3p"),
        "cap+": _make_evidence("cap", "3p", has_terminal_telomere=True),
    }

    resolver = RdnaCassetteResolver(graph, ledger, annotations=annotations)
    res = resolver.resolve_rdna_cassette(chrom="II", side="5p", entry_port=p_stem)

    assert res.end_state_record.resolved is True
    assert res.state == EndState.END_SHARED_EXEMPT
    assert res.quirk_record.quirk_type == IntervalQuirk.RDNA_CASSETTE

    # Test TSV emission via helper methods and standalone emitters
    end_tsv = res.emit_end_states_tsv()
    quirk_tsv = res.emit_quirks_tsv()

    assert "end_shared_exempt" in end_tsv
    assert "rdna_cassette" in quirk_tsv

    parsed_ends = parse_end_states_tsv(end_tsv)
    parsed_quirks = parse_quirks_tsv(StringIO(quirk_tsv))

    assert len(parsed_ends) == 1
    assert parsed_ends[0].state == EndState.END_SHARED_EXEMPT
    assert parsed_ends[0].shared_port == "cap+"
    assert parsed_ends[0].resolved is True

    assert len(parsed_quirks) == 1
    assert parsed_quirks[0].quirk_type == IntervalQuirk.RDNA_CASSETTE
    assert parsed_quirks[0].left_anchor_port == "stem+"
    assert parsed_quirks[0].right_anchor_port == "cap+"


def test_tuple_unpacking_support():
    """Verify RdnaResolution supports direct tuple unpacking into (end_state_record, quirk_record)."""
    graph = BidirectedGraph()
    graph.add_segment("stem", length=20000, depth=60.0)
    graph.add_segment("rdna", length=10000, depth=200.0)
    graph.add_segment("cap", length=5000, depth=30.0)

    p_stem = SegmentPort("stem", Strand.FORWARD)
    p_rdna = SegmentPort("rdna", Strand.FORWARD)
    p_cap = SegmentPort("cap", Strand.FORWARD)

    graph.add_edge(p_stem, p_rdna, cigar="100M")
    graph.add_edge(p_rdna, p_cap, cigar="100M")

    annotations = {
        "rdna": NodeAnnotations(
            sequence_id="rdna",
            barrnap_hits=[BarrnapHit("rdna", 100, 1800, "+", "18S_rRNA", 1e-40)],
        ),
    }

    ledger = {
        "stem+": _make_evidence("stem", "3p"),
        "cap+": _make_evidence("cap", "3p", has_terminal_telomere=True),
    }

    resolver = RdnaCassetteResolver(graph, ledger, annotations=annotations)
    end_rec, quirk_rec = resolver.resolve_rdna_cassette(chrom="II", side="5p", entry_port=p_stem)

    assert end_rec.state == EndState.END_SHARED_EXEMPT
    assert end_rec.resolved is True
    assert quirk_rec.quirk_type == IntervalQuirk.RDNA_CASSETTE


@pytest.mark.bio
@pytest.mark.benchmark
def test_benchmark_df5120_rdna_traversal_chr_ii():
    """Verify traversal from Chr II core anchor utg000013l through rDNA tangle to distal cap utg000053l-."""
    if not BENCHMARK_DIR.exists():
        pytest.skip(f"Benchmark DF5120 dataset not present in {BENCHMARK_DIR}")

    gfa_candidates = [
        BENCHMARK_DIR / "asm.p_utg.gfa.gz",
        BENCHMARK_DIR / "asm.p_utg.gfa",
    ]
    gfa_path = next((p for p in gfa_candidates if p.exists()), None)
    ledger_path = BENCHMARK_DIR / "df5120_end_evidence_ledger.tsv"
    barrnap_path = BENCHMARK_DIR / "barrnap_rdna.gff"

    if gfa_path is None or not ledger_path.exists() or not barrnap_path.exists():
        pytest.skip(f"Benchmark files missing in {BENCHMARK_DIR}")

    graph = load_gfa(gfa_path)
    barrnap = load_barrnap(barrnap_path)

    # Convert barrnap into NodeAnnotations
    annotations: dict[str, NodeAnnotations] = {}
    for seg_name, hits in barrnap.items():
        annotations[seg_name] = NodeAnnotations(sequence_id=seg_name, barrnap_hits=hits)

    # Load ledger
    ledger: dict[str, EndEvidence] = {}
    with open(ledger_path, encoding="utf-8") as f:
        _ = next(f)
        for line in f:
            row = line.strip().split("\t")
            if not row or len(row) < 17:
                continue
            seg_name = row[0]
            end = row[1]
            port_sym = row[2]
            length = int(row[3])
            depth = float(row[4])
            degree = int(row[5])
            is_tip = row[6] == "True"
            out_neigh = tuple(row[7].split(",")) if row[7] else ()
            has_tel = row[12] == "True"
            rdna_count = int(row[13])
            has_rdna = row[15] == "True"
            usco_count = int(row[16]) if len(row) > 16 and row[16] else 0
            strand = Strand.MINUS if end == "5p" else Strand.PLUS
            port = SegmentPort(seg_name, strand)
            ev = EndEvidence(
                segment_name=seg_name,
                end=end,
                port=port,
                length=length,
                depth=depth,
                degree=degree,
                is_tip=is_tip,
                outgoing_neighbors=out_neigh,
                has_terminal_telomere=has_tel,
                rdna_hit_count=rdna_count,
                has_terminal_rdna=has_rdna,
                usco_hit_count=usco_count,
            )
            ledger[port_sym] = ev
            ledger[f"{seg_name}:{end}"] = ev

    # Connect distal cap utg000053l- to rDNA tangle if not present in primary unitig graph
    # (utg000053l- has terminal telomere repeats and connects to rDNA unitigs in string graph)
    p_cap = SegmentPort("utg000053l", Strand.REVERSE)
    if not graph.edges.get(p_cap):
        # In Hifiasm string graph, utg000020l- connects to distal cap utg000053l-
        graph.add_edge(SegmentPort("utg000020l", Strand.REVERSE), p_cap, cigar="100M")

    resolver = RdnaCassetteResolver(graph, ledger, annotations=annotations)

    # 1. Entry ports must include Chr II core anchor utg000013l-
    entries = resolver.find_rdna_entry_ports()
    p_anchor = SegmentPort("utg000013l", Strand.REVERSE)
    assert p_anchor in entries

    # 2. Traversal reaches distal telomeric cap utg000053l-
    path = resolver.reach_distal_telomere_cap(p_anchor)
    assert path is not None
    assert path[0] == p_anchor
    assert path[-1] == p_cap
    assert any(resolver.is_rdna_unitig(p) for p in path)

    # 3. Resolve cassette creates valid records
    res = resolver.resolve_rdna_cassette(chrom="II", side="5p", entry_port=p_anchor)
    assert res.state == EndState.END_SHARED_EXEMPT
    assert res.end_state_record.resolved is True
    assert res.end_state_record.shared_port == "utg000053l-"
    assert res.quirk_record.quirk_type == IntervalQuirk.RDNA_CASSETTE
    assert res.quirk_record.left_anchor_port == "utg000013l-"
    assert res.quirk_record.right_anchor_port == "utg000053l-"


def test_distal_cap_not_misidentified_as_entry_port():
    """Verify distal telomeric cap is not misidentified as entry port due to bidirected duality."""
    graph = BidirectedGraph()
    graph.add_segment("z_stem", length=20000, depth=60.0, sequence="A" * 20000)
    graph.add_segment("rdna_1", length=10000, depth=200.0, sequence="C" * 10000)
    # Alphabetically earlier name than z_stem
    graph.add_segment("a_cap", length=5000, depth=30.0, sequence="T" * 5000)

    p_stem = SegmentPort("z_stem", Strand.FORWARD)
    p_r1 = SegmentPort("rdna_1", Strand.FORWARD)
    p_cap = SegmentPort("a_cap", Strand.FORWARD)

    graph.add_edge(p_stem, p_r1, cigar="100M")
    graph.add_edge(p_r1, p_cap, cigar="100M")

    annotations = {
        "rdna_1": NodeAnnotations(
            sequence_id="rdna_1",
            barrnap_hits=[
                BarrnapHit("rdna_1", 100, 1800, "+", "18S_rRNA", 1e-40),
                BarrnapHit("rdna_1", 2500, 6000, "+", "28S_rRNA", 1e-60),
            ],
        ),
    }

    ledger = {
        "z_stem+": _make_evidence("z_stem", "3p"),
        "a_cap+": _make_evidence("a_cap", "3p", has_terminal_telomere=True),
        "a_cap-": _make_evidence("a_cap", "5p", has_terminal_telomere=False),
        "rdna_1+": _make_evidence("rdna_1", "3p", rdna_hit_count=3),
    }

    resolver = RdnaCassetteResolver(graph, ledger, annotations=annotations)

    # find_rdna_entry_ports must NOT include a_cap-
    entry_ports = resolver.find_rdna_entry_ports()
    assert entry_ports == [p_stem]
    assert resolver.identify_entry_anchor() == p_stem

    # Automatic resolution without explicit entry_port must succeed
    res = resolver.resolve_rdna_cassette()
    assert res.end_state_record.resolved is True
    assert res.end_state_record.shared_port == "a_cap+"
    assert res.quirk_record.left_anchor_port == "z_stem+"
    assert res.quirk_record.right_anchor_port == "a_cap+"


def test_distal_cap_wrong_strand_not_accepted():
    """Verify distal telomeric cap entered with non-telomeric exit strand is rejected."""
    graph = BidirectedGraph()
    graph.add_segment("stem", length=20000, depth=60.0)
    graph.add_segment("rdna", length=10000, depth=200.0)
    graph.add_segment("cap", length=5000, depth=30.0)

    p_stem = SegmentPort("stem", Strand.FORWARD)
    p_rdna = SegmentPort("rdna", Strand.FORWARD)
    # cap has telomere only at 5p (cap-). Edge connects to cap+ (3p), which is naked.
    p_cap_naked = SegmentPort("cap", Strand.FORWARD)

    graph.add_edge(p_stem, p_rdna, cigar="100M")
    graph.add_edge(p_rdna, p_cap_naked, cigar="100M")

    annotations = {
        "rdna": NodeAnnotations(
            sequence_id="rdna",
            barrnap_hits=[BarrnapHit("rdna", 100, 1800, "+", "18S_rRNA", 1e-40)],
        ),
    }

    ledger = {
        "stem+": _make_evidence("stem", "3p"),
        "cap+": _make_evidence("cap", "3p", has_terminal_telomere=False),
        "cap-": _make_evidence("cap", "5p", has_terminal_telomere=True),
        "rdna+": _make_evidence("rdna", "3p", rdna_hit_count=1),
    }

    resolver = RdnaCassetteResolver(graph, ledger, annotations=annotations)
    assert resolver.is_distal_telomere_cap(p_cap_naked) is False
    assert resolver.is_distal_telomere_cap(SegmentPort("cap", Strand.REVERSE)) is True

    # Path to cap+ cannot be accepted as distal telomere cap
    assert resolver.reach_distal_telomere_cap(p_stem) is None


def test_internal_inversion_switch_rejected():
    """Verify unitigs with internal inverted rRNA orders (28S before 18S) or mixed strands are rejected."""
    graph = BidirectedGraph()
    graph.add_segment("stem", length=20000, depth=60.0)
    graph.add_segment("rdna_inverted", length=10000, depth=200.0)
    graph.add_segment("rdna_chimera", length=10000, depth=200.0)

    p_stem = SegmentPort("stem", Strand.FORWARD)
    p_inv = SegmentPort("rdna_inverted", Strand.FORWARD)
    p_chim = SegmentPort("rdna_chimera", Strand.FORWARD)

    graph.add_edge(p_stem, p_inv, cigar="100M")
    graph.add_edge(p_stem, p_chim, cigar="100M")

    annotations = {
        # Inverted: 28S at 100..2000, 18S at 3000..4800 on + strand
        "rdna_inverted": NodeAnnotations(
            sequence_id="rdna_inverted",
            barrnap_hits=[
                BarrnapHit("rdna_inverted", 100, 2000, "+", "28S_rRNA", 1e-50),
                BarrnapHit("rdna_inverted", 3000, 4800, "+", "18S_rRNA", 1e-50),
            ],
        ),
        # Chimera: 18S on +, 28S on -
        "rdna_chimera": NodeAnnotations(
            sequence_id="rdna_chimera",
            barrnap_hits=[
                BarrnapHit("rdna_chimera", 100, 1800, "+", "18S_rRNA", 1e-50),
                BarrnapHit("rdna_chimera", 2500, 6000, "-", "28S_rRNA", 1e-50),
            ],
        ),
    }

    ledger = {
        "stem+": _make_evidence("stem", "3p"),
        "rdna_inverted+": _make_evidence("rdna_inverted", "3p", rdna_hit_count=2),
        "rdna_chimera+": _make_evidence("rdna_chimera", "3p", rdna_hit_count=2),
    }

    resolver = RdnaCassetteResolver(graph, ledger, annotations=annotations)

    # Inverted unitig has no valid intrinsic orientation and cannot be traversed on either strand
    assert resolver.get_intrinsic_rdna_orientation("rdna_inverted") is None
    assert resolver.is_consistent_rdna_orientation(p_inv) is False
    assert resolver.is_consistent_rdna_orientation(SegmentPort("rdna_inverted", Strand.REVERSE)) is False

    # Chimeric unitig is likewise rejected
    assert resolver.get_intrinsic_rdna_orientation("rdna_chimera") is None
    assert resolver.is_consistent_rdna_orientation(p_chim) is False
    assert resolver.is_consistent_rdna_orientation(SegmentPort("rdna_chimera", Strand.REVERSE)) is False


def test_entry_port_ranking_prefers_usco_core_stem():
    """Verify entry anchor selection prioritizes chromosome core stem with USCO genes over debris."""
    graph = BidirectedGraph()
    graph.add_segment("a_debris", length=500, depth=20.0)
    graph.add_segment("z_stem", length=50000, depth=60.0)
    graph.add_segment("rdna", length=10000, depth=200.0)

    p_debris = SegmentPort("a_debris", Strand.FORWARD)
    p_stem = SegmentPort("z_stem", Strand.FORWARD)
    p_rdna = SegmentPort("rdna", Strand.FORWARD)

    graph.add_edge(p_debris, p_rdna, cigar="100M")
    graph.add_edge(p_stem, p_rdna, cigar="100M")

    annotations = {
        "rdna": NodeAnnotations(
            sequence_id="rdna",
            barrnap_hits=[BarrnapHit("rdna", 100, 1800, "+", "18S_rRNA", 1e-40)],
        ),
    }

    ledger = {
        "a_debris+": _make_evidence("a_debris", "3p", length=500, usco_hit_count=0),
        "z_stem+": _make_evidence("z_stem", "3p", length=50000, usco_hit_count=15),
        "rdna+": _make_evidence("rdna", "3p", rdna_hit_count=1),
    }

    resolver = RdnaCassetteResolver(graph, ledger, annotations=annotations)

    # z_stem must be ranked first despite a_debris coming first alphabetically
    entry_ports = resolver.find_rdna_entry_ports()
    assert entry_ports[0] == p_stem
    assert resolver.identify_entry_anchor() == p_stem


def test_direct_stem_to_cap_edge_without_rdna_not_accepted():
    """Verify direct connection from stem to cap without rDNA is not accepted as rDNA cassette."""
    graph = BidirectedGraph()
    graph.add_segment("stem", length=20000, depth=60.0)
    graph.add_segment("cap", length=5000, depth=30.0)

    p_stem = SegmentPort("stem", Strand.FORWARD)
    p_cap = SegmentPort("cap", Strand.FORWARD)

    # Direct stem -> cap edge (no rDNA)
    graph.add_edge(p_stem, p_cap, cigar="100M")

    ledger = {
        "stem+": _make_evidence("stem", "3p"),
        "cap+": _make_evidence("cap", "3p", has_terminal_telomere=True),
    }

    resolver = RdnaCassetteResolver(graph, ledger)
    assert resolver.reach_distal_telomere_cap(p_stem) is None
    assert resolver.traverse_rdna_cassette(p_stem) == []


def test_compute_path_length_matches_assembled_sequence_with_indels():
    """Verify compute_path_length accurately reflects assembled sequence length with CIGAR indels."""
    graph = BidirectedGraph()
    # seg1: 1000 bp, seg2: 500 bp. CIGAR has 2 bp insertion (query consumes 100 bp, target consumes 98 bp)
    seq1 = "A" * 1000
    seq2 = "C" * 500
    graph.add_segment("s1", length=1000, sequence=seq1)
    graph.add_segment("s2", length=500, sequence=seq2)

    p1 = SegmentPort("s1", Strand.FORWARD)
    p2 = SegmentPort("s2", Strand.FORWARD)

    graph.add_edge(p1, p2, cigar="50M2I48M")

    resolver = RdnaCassetteResolver(graph, {})
    path = [p1, p2]

    assembled = resolver.assemble_contiguous_path(path)
    assert assembled is not None
    assert len(assembled) == resolver.compute_path_length(path)


@pytest.mark.bio
@pytest.mark.benchmark
def test_benchmark_df5120_rdna_chains_eukaryotic_and_rejects_bacterial_contaminants():
    """Verify RdnaCassetteResolver chains utg000013l- -> utg000020l- -> utg000082l- and rejects bacterial contigs."""
    if not BENCHMARK_DIR.exists():
        pytest.skip(f"Benchmark DF5120 dataset not present in {BENCHMARK_DIR}")

    gfa_candidates = [
        BENCHMARK_DIR / "asm.p_utg.gfa.gz",
        BENCHMARK_DIR / "asm.p_utg.gfa",
    ]
    gfa_path = next((p for p in gfa_candidates if p.exists()), None)
    if gfa_path is None:
        pytest.skip(f"GFA file missing in {BENCHMARK_DIR}")

    tiara_path = BENCHMARK_DIR / "tiara_unitigs.txt"
    tidk_path = BENCHMARK_DIR / "tidk_telomeres_telomeric_repeat_windows.tsv"
    barrnap_path = BENCHMARK_DIR / "barrnap_rdna.gff"
    compleasm_path = BENCHMARK_DIR / "compleasm_unitigs_full_table.tsv"

    if not (tiara_path.exists() and tidk_path.exists() and barrnap_path.exists() and compleasm_path.exists()):
        pytest.skip(f"Benchmark annotation files missing in {BENCHMARK_DIR}")

    from gfapdefinisher.annotations import aggregate_node_annotations
    from gfapdefinisher.ledger import build_end_evidence_ledger

    graph = load_gfa(gfa_path)
    annotations = aggregate_node_annotations(
        tiara=tiara_path,
        tidk=tidk_path,
        barrnap=barrnap_path,
        compleasm=compleasm_path,
        sequence_ids=graph.segments.keys(),
    )
    ledger = build_end_evidence_ledger(graph=graph, annotations=annotations)
    ledger_map = {ev.port.to_symbol(): ev for ev in ledger} | {f"{ev.segment_name}:{ev.end}": ev for ev in ledger}

    resolver = RdnaCassetteResolver(graph, ledger_map, annotations=annotations)

    # 1. Contaminant detection checks
    assert resolver._is_contaminant("utg000005l") is True
    assert resolver._is_contaminant("utg000008l") is True
    assert resolver._is_contaminant("utg000013l") is False
    assert resolver._is_contaminant("utg000020l") is False
    assert resolver._is_contaminant("utg000082l") is False

    # 2. Entry port detection excludes bacterial contigs
    entries = resolver.find_rdna_entry_ports()
    p_anchor = SegmentPort("utg000013l", Strand.REVERSE)
    assert p_anchor in entries
    entry_names = {p.name for p in entries}
    banned_bacterial = {"utg000005l", "utg000008l", "utg000075l", "utg000156l", "utg000192l"}
    for b in banned_bacterial:
        assert b not in entry_names

    # 3. Traversal chains eukaryotic rDNA copies utg000013l- -> utg000020l- -> utg000082l- -> utg000214l+ -> utg000212l+ -> utg000172l+ -> utg000207l-
    chain = resolver.traverse_rdna_cassette(p_anchor, include_anchor=True)
    chain_symbols = [p.to_symbol() for p in chain]
    assert chain_symbols == [
        "utg000013l-",
        "utg000020l-",
        "utg000082l-",
        "utg000214l+",
        "utg000212l+",
        "utg000172l+",
        "utg000207l-",
    ]
    assert not any(p.name in banned_bacterial for p in chain)


def test_resolve_entry_stem_heuristic_orphan_cap_matching_and_rejection():
    """Verify resolve_entry_stem populates distal_cap_port and capped_via_heuristic=True on matching orphan cap,

    and rejects orphan caps with opposite rDNA orientation or contaminant classification.
    """
    # 1. Matching orphan cap test
    graph = BidirectedGraph()
    graph.add_segment("stem", length=20000, depth=60.0)
    graph.add_segment("rdna_tangle_end", length=10000, depth=200.0)
    graph.add_segment("orphan_cap", length=5000, depth=30.0)

    p_stem = SegmentPort("stem", Strand.FORWARD)
    p_tangle = SegmentPort("rdna_tangle_end", Strand.FORWARD)
    graph.add_edge(p_stem, p_tangle, cigar="100M")
    # orphan_cap is completely disconnected (degree 0 at both ports)

    annotations = {
        "rdna_tangle_end": NodeAnnotations(
            sequence_id="rdna_tangle_end",
            barrnap_hits=[
                BarrnapHit("rdna_tangle_end", 100, 1800, "+", "18S_rRNA", 1e-40),
                BarrnapHit("rdna_tangle_end", 2500, 6000, "+", "28S_rRNA", 1e-60),
            ],
        ),
        "orphan_cap": NodeAnnotations(
            sequence_id="orphan_cap",
            barrnap_hits=[
                BarrnapHit("orphan_cap", 100, 1800, "+", "18S_rRNA", 1e-40),
                BarrnapHit("orphan_cap", 2500, 4500, "+", "28S_rRNA", 1e-60),
            ],
        ),
    }

    ledger = {
        "stem+": _make_evidence("stem", "3p"),
        "rdna_tangle_end+": _make_evidence("rdna_tangle_end", "3p", rdna_hit_count=2),
        "orphan_cap+": _make_evidence("orphan_cap", "3p", has_terminal_telomere=True, rdna_hit_count=2),
        "orphan_cap-": _make_evidence("orphan_cap", "5p", has_terminal_telomere=False, rdna_hit_count=2),
    }

    resolver = RdnaCassetteResolver(graph, ledger, annotations=annotations)
    res = resolver.resolve_entry_stem(p_stem, chrom="II", side="5p")

    assert res.distal_cap_port == SegmentPort("orphan_cap", Strand.FORWARD)
    assert res.capped_via_heuristic is True
    assert res.reached_distal_cap is False
    assert res.quirk_record.agp_component_type == "U"
    assert res.quirk_record.agp_linkage_evidence == "rna_tangle"
    assert res.quirk_record.right_anchor_port == "orphan_cap+"
    assert res.end_state_record.resolved is False

    # 2. Rejection of opposite rDNA orientation
    graph_opp = BidirectedGraph()
    graph_opp.add_segment("stem", length=20000, depth=60.0)
    graph_opp.add_segment("rdna_tangle_end", length=10000, depth=200.0)
    graph_opp.add_segment("orphan_cap_opp", length=5000, depth=30.0)
    graph_opp.add_edge(p_stem, p_tangle, cigar="100M")

    annotations_opp = {
        "rdna_tangle_end": NodeAnnotations(
            sequence_id="rdna_tangle_end",
            barrnap_hits=[
                BarrnapHit("rdna_tangle_end", 100, 1800, "+", "18S_rRNA", 1e-40),
                BarrnapHit("rdna_tangle_end", 2500, 6000, "+", "28S_rRNA", 1e-60),
            ],
        ),
        # Opposite orientation: intrinsic minus strand (decreasing coordinates: 18S at 4500, 28S at 100)
        "orphan_cap_opp": NodeAnnotations(
            sequence_id="orphan_cap_opp",
            barrnap_hits=[
                BarrnapHit("orphan_cap_opp", 4000, 4800, "-", "18S_rRNA", 1e-40),
                BarrnapHit("orphan_cap_opp", 500, 2500, "-", "28S_rRNA", 1e-60),
            ],
        ),
    }

    ledger_opp = {
        "stem+": _make_evidence("stem", "3p"),
        "rdna_tangle_end+": _make_evidence("rdna_tangle_end", "3p", rdna_hit_count=2),
        "orphan_cap_opp+": _make_evidence("orphan_cap_opp", "3p", has_terminal_telomere=True, rdna_hit_count=2),
        "orphan_cap_opp-": _make_evidence("orphan_cap_opp", "5p", has_terminal_telomere=False, rdna_hit_count=2),
    }

    resolver_opp = RdnaCassetteResolver(graph_opp, ledger_opp, annotations=annotations_opp)
    res_opp = resolver_opp.resolve_entry_stem(p_stem, chrom="II", side="5p")

    assert res_opp.distal_cap_port is None
    assert res_opp.capped_via_heuristic is False
    assert res_opp.quirk_record.right_anchor_port is None

    # 3. Rejection of contaminant classification
    graph_contam = BidirectedGraph()
    graph_contam.add_segment("stem", length=20000, depth=60.0)
    graph_contam.add_segment("rdna_tangle_end", length=10000, depth=200.0)
    graph_contam.add_segment("orphan_cap_contam", length=5000, depth=30.0)
    graph_contam.add_edge(p_stem, p_tangle, cigar="100M")

    annotations_contam = {
        "rdna_tangle_end": NodeAnnotations(
            sequence_id="rdna_tangle_end",
            barrnap_hits=[
                BarrnapHit("rdna_tangle_end", 100, 1800, "+", "18S_rRNA", 1e-40),
                BarrnapHit("rdna_tangle_end", 2500, 6000, "+", "28S_rRNA", 1e-60),
            ],
        ),
        "orphan_cap_contam": NodeAnnotations(
            sequence_id="orphan_cap_contam",
            barrnap_hits=[
                BarrnapHit("orphan_cap_contam", 100, 1800, "+", "18S_rRNA", 1e-40),
                BarrnapHit("orphan_cap_contam", 2500, 4500, "+", "28S_rRNA", 1e-60),
            ],
            tiara="bacteria",
        ),
    }

    ledger_contam = {
        "stem+": _make_evidence("stem", "3p"),
        "rdna_tangle_end+": _make_evidence("rdna_tangle_end", "3p", rdna_hit_count=2),
        "orphan_cap_contam+": _make_evidence("orphan_cap_contam", "3p", has_terminal_telomere=True, rdna_hit_count=2),
        "orphan_cap_contam-": _make_evidence("orphan_cap_contam", "5p", has_terminal_telomere=False, rdna_hit_count=2),
    }

    resolver_contam = RdnaCassetteResolver(graph_contam, ledger_contam, annotations=annotations_contam)
    res_contam = resolver_contam.resolve_entry_stem(p_stem, chrom="II", side="5p")

    assert res_contam.distal_cap_port is None
    assert res_contam.capped_via_heuristic is False
    assert res_contam.quirk_record.right_anchor_port is None


def test_resolve_entry_stem_concatenates_all_rdna_unitigs_and_records_excluded_caps():
    """Verify resolve_entry_stem greedily chains disconnected rDNA unitigs and records excluded orphan caps."""
    graph = BidirectedGraph()
    graph.add_segment("stem", length=20000, depth=60.0)
    graph.add_segment("rdna_1", length=10000, depth=100.0)
    graph.add_segment("rdna_2", length=8000, depth=80.0)
    graph.add_segment("rdna_3", length=6000, depth=40.0)
    graph.add_segment("cap_high", length=5000, depth=30.0)
    graph.add_segment("cap_low", length=4000, depth=5.0)

    p_stem = SegmentPort("stem", Strand.FORWARD)
    p_rdna_1 = SegmentPort("rdna_1", Strand.FORWARD)
    graph.add_edge(p_stem, p_rdna_1, cigar="100M")

    # rdna_2 and rdna_3 are disconnected from rdna_1, but rdna_2 connects to rdna_3
    p_rdna_2 = SegmentPort("rdna_2", Strand.FORWARD)
    p_rdna_3 = SegmentPort("rdna_3", Strand.FORWARD)
    graph.add_edge(p_rdna_2, p_rdna_3, cigar="100M")

    annotations = {
        "rdna_1": NodeAnnotations(
            sequence_id="rdna_1",
            barrnap_hits=[
                BarrnapHit("rdna_1", 100, 1800, "+", "18S_rRNA", 1e-40),
                BarrnapHit("rdna_1", 2500, 6000, "+", "28S_rRNA", 1e-60),
            ],
        ),
        "rdna_2": NodeAnnotations(
            sequence_id="rdna_2",
            barrnap_hits=[
                BarrnapHit("rdna_2", 100, 1800, "+", "18S_rRNA", 1e-40),
                BarrnapHit("rdna_2", 2500, 6000, "+", "28S_rRNA", 1e-60),
            ],
        ),
        "rdna_3": NodeAnnotations(
            sequence_id="rdna_3",
            barrnap_hits=[
                BarrnapHit("rdna_3", 100, 1800, "+", "18S_rRNA", 1e-40),
                BarrnapHit("rdna_3", 2500, 5000, "+", "28S_rRNA", 1e-60),
            ],
        ),
        "cap_high": NodeAnnotations(
            sequence_id="cap_high",
            barrnap_hits=[
                BarrnapHit("cap_high", 100, 1800, "+", "18S_rRNA", 1e-40),
                BarrnapHit("cap_high", 2500, 4500, "+", "28S_rRNA", 1e-60),
            ],
        ),
        "cap_low": NodeAnnotations(
            sequence_id="cap_low",
            barrnap_hits=[
                BarrnapHit("cap_low", 100, 1800, "+", "18S_rRNA", 1e-40),
            ],
        ),
    }

    ledger = {
        "stem+": _make_evidence("stem", "3p"),
        "rdna_1+": _make_evidence("rdna_1", "3p", rdna_hit_count=2),
        "rdna_2+": _make_evidence("rdna_2", "3p", rdna_hit_count=2),
        "rdna_3+": _make_evidence("rdna_3", "3p", rdna_hit_count=2),
        "cap_high+": _make_evidence("cap_high", "3p", has_terminal_telomere=True, rdna_hit_count=2),
        "cap_high-": _make_evidence("cap_high", "5p", has_terminal_telomere=False, rdna_hit_count=2),
        "cap_low+": _make_evidence("cap_low", "3p", has_terminal_telomere=True, rdna_hit_count=1),
        "cap_low-": _make_evidence("cap_low", "5p", has_terminal_telomere=False, rdna_hit_count=1),
    }

    resolver = RdnaCassetteResolver(graph, ledger, annotations=annotations)
    res = resolver.resolve_entry_stem(p_stem, chrom="II", side="5p")

    # Path must chain stem, rdna_1, then jump to rdna_2, then step to rdna_3
    assert [p.to_symbol() for p in res.path] == ["stem+", "rdna_1+", "rdna_2+", "rdna_3+"]
    assert res.unlinked_join_indices == (2,)  # jump at rdna_2
    assert res.distal_cap_port == SegmentPort("cap_high", Strand.FORWARD)
    assert res.excluded_unitigs == ("cap_low",)
    assert res.quirk_record.unitig_ids == ("stem", "rdna_1", "rdna_2", "rdna_3", "cap_high")


@pytest.mark.bio
@pytest.mark.benchmark
def test_benchmark_df5120_full_rdna_concatenation_and_cap_exclusion():
    """Verify on reference benchmark DF5120 that RdnaCassetteResolver chains all 16 rDNA unitigs and caps with utg000053l-."""
    if not BENCHMARK_DIR.exists():
        pytest.skip(f"Benchmark DF5120 dataset not present in {BENCHMARK_DIR}")

    gfa_candidates = [
        BENCHMARK_DIR / "asm.p_utg.gfa.gz",
        BENCHMARK_DIR / "asm.p_utg.gfa",
    ]
    gfa_path = next((p for p in gfa_candidates if p.exists()), None)
    if gfa_path is None:
        pytest.skip(f"GFA file missing in {BENCHMARK_DIR}")

    from gfapdefinisher.annotations import aggregate_node_annotations
    from gfapdefinisher.ledger import build_end_evidence_ledger

    graph = load_gfa(gfa_path)
    annotations = aggregate_node_annotations(
        tiara=BENCHMARK_DIR / "tiara_unitigs.txt",
        tidk=BENCHMARK_DIR / "tidk_telomeres_telomeric_repeat_windows.tsv",
        barrnap=BENCHMARK_DIR / "barrnap_rdna.gff",
        compleasm=BENCHMARK_DIR / "compleasm_unitigs_full_table.tsv",
        sequence_ids=graph.segments.keys(),
    )
    ledger = build_end_evidence_ledger(graph=graph, annotations=annotations)
    ledger_map = {ev.port.to_symbol(): ev for ev in ledger} | {f"{ev.segment_name}:{ev.end}": ev for ev in ledger}

    resolver = RdnaCassetteResolver(graph, ledger_map, annotations=annotations)
    p_anchor = SegmentPort("utg000013l", Strand.REVERSE)
    res = resolver.resolve_entry_stem(p_anchor, chrom="II", side="5p")

    # Chains core stem + 15 rDNA tangle unitigs
    assert len(res.path) == 16
    assert res.path[0] == p_anchor
    assert res.distal_cap_port == SegmentPort("utg000053l", Strand.REVERSE)
    assert res.excluded_unitigs == ("utg000098l",)
    assert len(res.quirk_record.unitig_ids) == 17
    assert res.quirk_record.right_anchor_port == "utg000053l-"
    assert res.quirk_record.agp_linkage_evidence == "rna_tangle"


def test_is_rdna_unitig_rejects_5s_only():
    """Verify is_rdna_unitig rejects 5S rRNA and prokaryotic rRNA (16S, 23S) hits, accepting only 45S."""
    hit_5s = BarrnapHit(
        sequence_id="utg_5s",
        start=1,
        end=120,
        strand="+",
        feature_name="5S_rRNA",
        evalue=1e-10,
    )
    ann_5s = {"utg_5s": NodeAnnotations(sequence_id="utg_5s", barrnap_hits=(hit_5s,))}
    assert is_rdna_unitig("utg_5s", annotations=ann_5s) is False

    # Ledger evidence with only 5S
    ev_5s = _make_evidence("utg_5s_ev", "3p", rdna_hit_count=1, rdna_features=("5S_rRNA",))
    assert is_rdna_unitig("utg_5s_ev", ledger={"utg_5s_ev+": ev_5s}) is False

    # Prokaryotic rRNA hits
    hit_16s = BarrnapHit(
        sequence_id="utg_16s",
        start=1,
        end=1500,
        strand="+",
        feature_name="16S_rRNA",
        evalue=1e-50,
    )
    ann_16s = {"utg_16s": NodeAnnotations(sequence_id="utg_16s", barrnap_hits=(hit_16s,))}
    assert is_rdna_unitig("utg_16s", annotations=ann_16s) is False

    hit_23s = BarrnapHit(
        sequence_id="utg_23s",
        start=1,
        end=2900,
        strand="+",
        feature_name="23S_rRNA",
        evalue=1e-50,
    )
    ann_23s = {"utg_23s": NodeAnnotations(sequence_id="utg_23s", barrnap_hits=(hit_23s,))}
    assert is_rdna_unitig("utg_23s", annotations=ann_23s) is False

    # Valid eukaryotic 45S rRNA hits
    hit_18s = BarrnapHit(
        sequence_id="utg_18s",
        start=1,
        end=1800,
        strand="+",
        feature_name="18S_rRNA",
        evalue=1e-50,
    )
    ann_18s = {"utg_18s": NodeAnnotations(sequence_id="utg_18s", barrnap_hits=(hit_18s,))}
    assert is_rdna_unitig("utg_18s", annotations=ann_18s) is True

    hit_28s = BarrnapHit(
        sequence_id="utg_28s",
        start=1,
        end=3500,
        strand="+",
        feature_name="28S_rRNA",
        evalue=1e-50,
    )
    ann_28s = {"utg_28s": NodeAnnotations(sequence_id="utg_28s", barrnap_hits=(hit_28s,))}
    assert is_rdna_unitig("utg_28s", annotations=ann_28s) is True


def test_intercalary_45s_rdna_cluster_bridging_core_stems():
    """Verify intercalary 45S rDNA cluster bridges two core chromosome stems without terminal telomeres."""
    graph = BidirectedGraph()
    graph.add_segment("stem_A", length=50000, depth=60.0, sequence="A" * 50000)
    graph.add_segment("rdna_int", length=15000, depth=180.0, sequence="C" * 15000)
    graph.add_segment("stem_B", length=60000, depth=60.0, sequence="G" * 60000)

    p_stem_A = SegmentPort("stem_A", Strand.FORWARD)
    p_rdna_int = SegmentPort("rdna_int", Strand.FORWARD)
    p_stem_B = SegmentPort("stem_B", Strand.FORWARD)

    # stem_A+ -> rdna_int+ -> stem_B+
    graph.add_edge(p_stem_A, p_rdna_int, cigar="100M")
    graph.add_edge(p_rdna_int, p_stem_B, cigar="100M")

    # Stems have USCO hits > 0 (core chromosomes), neither has terminal telomere
    ledger = {
        "stem_A+": _make_evidence("stem_A", "3p", usco_hit_count=5),
        "stem_A-": _make_evidence("stem_A", "5p", usco_hit_count=5),
        "rdna_int+": _make_evidence("rdna_int", "3p", rdna_hit_count=2, rdna_features=("18S_rRNA", "28S_rRNA")),
        "rdna_int-": _make_evidence("rdna_int", "5p", rdna_hit_count=2, rdna_features=("18S_rRNA", "28S_rRNA")),
        "stem_B+": _make_evidence("stem_B", "3p", usco_hit_count=6),
        "stem_B-": _make_evidence("stem_B", "5p", usco_hit_count=6),
    }

    h18 = BarrnapHit(sequence_id="rdna_int", start=1000, end=2800, strand="+", feature_name="18S_rRNA", evalue=1e-50)
    h28 = BarrnapHit(sequence_id="rdna_int", start=4000, end=7500, strand="+", feature_name="28S_rRNA", evalue=1e-50)
    ann = {"rdna_int": NodeAnnotations(sequence_id="rdna_int", barrnap_hits=(h18, h28))}

    resolver = RdnaCassetteResolver(graph, ledger, annotations=ann)
    res = resolver.resolve_entry_stem(p_stem_A, chrom="chrI", side="3p")

    assert res.is_intercalary is True
    assert res.reached_distal_cap is False
    assert res.exit_stem_port == p_stem_B
    assert [p.name for p in res.path] == ["stem_A", "rdna_int", "stem_B"]


def test_consumed_rdna_prevents_duplicate_cassette_traversal():
    """Verify consumed_rdna prevents duplicate traversal of the same 45S rDNA cluster from multiple entry ports."""
    graph = BidirectedGraph()
    graph.add_segment("stem_primary", length=80000, depth=60.0, sequence="A" * 80000)
    graph.add_segment("stem_secondary", length=20000, depth=60.0, sequence="T" * 20000)
    graph.add_segment("rdna_cluster", length=15000, depth=240.0, sequence="C" * 15000)
    graph.add_segment("tel_cap", length=5000, depth=30.0, sequence="G" * 5000)

    p_stem1 = SegmentPort("stem_primary", Strand.FORWARD)
    p_stem2 = SegmentPort("stem_secondary", Strand.FORWARD)
    p_rdna = SegmentPort("rdna_cluster", Strand.FORWARD)
    p_cap = SegmentPort("tel_cap", Strand.FORWARD)

    # Both stems connect into the same rDNA cluster
    graph.add_edge(p_stem1, p_rdna, cigar="100M")
    graph.add_edge(p_stem2, p_rdna, cigar="100M")
    graph.add_edge(p_rdna, p_cap, cigar="100M")

    ledger = {
        "stem_primary+": _make_evidence("stem_primary", "3p", usco_hit_count=10),
        "stem_primary-": _make_evidence("stem_primary", "5p", usco_hit_count=10),
        "stem_secondary+": _make_evidence("stem_secondary", "3p", usco_hit_count=2),
        "stem_secondary-": _make_evidence("stem_secondary", "5p", usco_hit_count=2),
        "rdna_cluster+": _make_evidence("rdna_cluster", "3p", rdna_hit_count=2, rdna_features=("18S_rRNA", "28S_rRNA")),
        "rdna_cluster-": _make_evidence("rdna_cluster", "5p", rdna_hit_count=2, rdna_features=("18S_rRNA", "28S_rRNA")),
        "tel_cap+": _make_evidence(
            "tel_cap",
            "3p",
            has_terminal_telomere=True,
            rdna_hit_count=1,
            rdna_features=("28S_rRNA",),
        ),
        "tel_cap-": _make_evidence(
            "tel_cap",
            "5p",
            has_terminal_telomere=False,
            rdna_hit_count=1,
            rdna_features=("28S_rRNA",),
        ),
    }

    h18 = BarrnapHit(
        sequence_id="rdna_cluster", start=1000, end=2800, strand="+", feature_name="18S_rRNA", evalue=1e-50
    )
    h28 = BarrnapHit(
        sequence_id="rdna_cluster", start=4000, end=7500, strand="+", feature_name="28S_rRNA", evalue=1e-50
    )
    ann = {"rdna_cluster": NodeAnnotations(sequence_id="rdna_cluster", barrnap_hits=(h18, h28))}

    resolver = RdnaCassetteResolver(graph, ledger, annotations=ann)
    consumed_rdna: set[str] = set()

    # Primary stem resolves the rDNA cluster
    res1 = resolver.resolve_entry_stem(p_stem1, chrom="chr1", side="3p", consumed_rdna=consumed_rdna)
    assert res1.reached_distal_cap or res1.capped_via_heuristic
    assert res1.distal_cap_port == p_cap
    for p in res1.path[1:]:
        consumed_rdna.add(p.name)
    if res1.distal_cap_port:
        consumed_rdna.add(res1.distal_cap_port.name)
    consumed_rdna.update(res1.excluded_unitigs)

    # Secondary stem entering the consumed cluster must return empty/unresolved, preventing duplicate cassette
    res2 = resolver.resolve_entry_stem(p_stem2, chrom="chr2", side="3p", consumed_rdna=consumed_rdna)
    assert len(res2.path) <= 1
    assert res2.distal_cap_port is None
    assert not res2.is_resolved
