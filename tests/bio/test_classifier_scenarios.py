"""Biological scenario verification tests for 3-axis classification and PDE junctions.

Falsification test suite validating all 12 authentic diminution boundaries
in nematode benchmark DF5120 (BIO-001, BIO-002, BIO-003, BIO-006):
  Megabase chromosome arm ends:
    1. utg000001l:3p (utg000001l+) -> somatic: utg000142l-, germline: utg000101l+
    2. utg000001l:5p (utg000001l-) -> somatic: utg000076l+, germline: utg000039l+
    3. utg000002l:3p (utg000002l+) -> somatic: utg000194l+, germline: utg000115l+
    4. utg000002l:5p (utg000002l-) -> somatic: utg000093l-, germline: utg000103l+
    5. utg000009l:5p (utg000009l-) -> somatic: utg000211l-, germline: utg000069l+
    6. utg000013l:3p (utg000013l+) -> somatic: utg000167l+, germline: utg000033l-
    7. utg000027l:5p (utg000027l-) -> somatic: utg000154l-, germline: utg000189l-
  Chr X internal breaks:
    8. utg000016l:5p (utg000016l-) -> somatic: utg000088l+, germline: utg000096l+
    9. utg000019l:3p (utg000019l+) -> somatic: utg000200l-, germline: utg000071l+
    10. utg000019l:5p (utg000019l-) -> somatic: utg000129l-, germline: utg000145l+
    11. utg000021l:5p (utg000021l-) -> somatic: utg000051l+, germline: utg000071l-
    12. utg000027l:3p (utg000027l+) -> somatic: utg000209l+, germline: utg000145l-
and rejecting negative controls utg000020l+ (rDNA) and utg000065l- (bubble).
"""

from pathlib import Path

import pytest

from gfapdefinisher.classifier import PdeJunctionCall, classify_pde_junction
from gfapdefinisher.graph import BidirectedGraph, load_gfa
from gfapdefinisher.ledger import EndEvidence
from gfapdefinisher.models import SegmentPort, Strand

BENCHMARK_DIR = Path("data/benchmark/DF5120")


def _load_benchmark_ledger(ledger_tsv_path: Path) -> dict[str, EndEvidence]:
    """Parse end evidence records from benchmark TSV ledger."""
    ledger: dict[str, EndEvidence] = {}
    with open(ledger_tsv_path, encoding="utf-8") as f:
        _ = next(f)  # Header
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
            tel_canon = int(row[8])
            tel_opp = int(row[9])
            has_tel = row[12] == "True"
            rdna_count = int(row[13])
            has_rdna = row[15] == "True"
            usco_count = int(row[16])
            tiara = row[18] if len(row) > 18 else ""

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
                telomere_canonical_count=tel_canon,
                telomere_opposite_count=tel_opp,
                has_terminal_telomere=has_tel,
                rdna_hit_count=rdna_count,
                has_terminal_rdna=has_rdna,
                usco_hit_count=usco_count,
                tiara_class=tiara,
            )
            ledger[port_sym] = ev
            ledger[f"{seg_name}:{end}"] = ev
    return ledger


@pytest.mark.bio
def test_synthetic_four_benchmark_pde_junction_scenarios():
    """Verify synthetic topology replicating all 4 benchmark diminution junctions."""
    benchmarks = [
        {
            "stem": ("utg000009l", Strand.MINUS, 62.0, "5p"),
            "som": ("utg000211l", Strand.MINUS, 27.0, True),
            "grs": ("utg000069l", Strand.FORWARD, 6.0, False),
            "expected_som": SegmentPort("utg000211l", Strand.MINUS),
            "expected_grs": SegmentPort("utg000069l", Strand.FORWARD),
        },
        {
            "stem": ("utg000016l", Strand.MINUS, 63.0, "5p"),
            "som": ("utg000088l", Strand.FORWARD, 45.0, True),
            "grs": ("utg000096l", Strand.FORWARD, 7.0, False),
            "expected_som": SegmentPort("utg000088l", Strand.FORWARD),
            "expected_grs": SegmentPort("utg000096l", Strand.FORWARD),
        },
        {
            "stem": ("utg000019l", Strand.FORWARD, 66.0, "3p"),
            "som": ("utg000200l", Strand.MINUS, 39.0, True),
            "grs": ("utg000071l", Strand.FORWARD, 7.0, False),
            "expected_som": SegmentPort("utg000200l", Strand.MINUS),
            "expected_grs": SegmentPort("utg000071l", Strand.FORWARD),
        },
        {
            "stem": ("utg000021l", Strand.MINUS, 64.0, "5p"),
            "som": ("utg000051l", Strand.FORWARD, 51.0, True),
            "grs": ("utg000071l", Strand.MINUS, 7.0, False),
            "expected_som": SegmentPort("utg000051l", Strand.FORWARD),
            "expected_grs": SegmentPort("utg000071l", Strand.MINUS),
        },
    ]

    for b in benchmarks:
        stem_name, stem_strand, stem_depth, stem_end = b["stem"]
        som_name, som_strand, som_depth, som_tel = b["som"]
        grs_name, grs_strand, grs_depth, grs_tel = b["grs"]

        graph = BidirectedGraph()
        graph.add_segment(stem_name, length=100000, depth=stem_depth)
        graph.add_segment(som_name, length=25000, depth=som_depth)
        graph.add_segment(grs_name, length=50000, depth=grs_depth)

        p_stem = SegmentPort(stem_name, stem_strand)
        p_som = SegmentPort(som_name, som_strand)
        p_grs = SegmentPort(grs_name, grs_strand)

        graph.add_edge(p_stem, p_som, cigar="100M")
        graph.add_edge(p_stem, p_grs, cigar="100M")

        ev_stem = EndEvidence(
            segment_name=stem_name,
            end=stem_end,
            port=p_stem,
            length=100000,
            depth=stem_depth,
            degree=2,
            is_tip=False,
            outgoing_neighbors=(p_som.to_symbol(), p_grs.to_symbol()),
            usco_hit_count=10,
            tiara_class="eukarya",
        )
        ev_som = EndEvidence(
            segment_name=som_name,
            end="3p" if som_strand == Strand.FORWARD else "5p",
            port=p_som,
            length=25000,
            depth=som_depth,
            degree=0,
            is_tip=True,
            outgoing_neighbors=(),
            has_terminal_telomere=som_tel,
            tiara_class="eukarya",
        )
        ev_grs = EndEvidence(
            segment_name=grs_name,
            end="3p" if grs_strand == Strand.FORWARD else "5p",
            port=p_grs,
            length=50000,
            depth=grs_depth,
            degree=0,
            is_tip=True,
            outgoing_neighbors=(),
            has_terminal_telomere=grs_tel,
            tiara_class="eukarya",
        )

        ledger = {
            p_stem.to_symbol(): ev_stem,
            p_som.to_symbol(): ev_som,
            p_grs.to_symbol(): ev_grs,
            f"{stem_name}:{stem_end}": ev_stem,
        }

        call: PdeJunctionCall = classify_pde_junction(p_stem, graph, ledger)
        assert call.is_pde_junction is True, f"Failed for {stem_name}:{stem_end}"
        assert call.port == p_stem
        assert call.somatic_port == b["expected_som"]
        assert call.germline_port == b["expected_grs"]
        assert 0.04 <= call.depth_ratio <= 0.33
        assert 0.04 <= call.core_depth_ratio <= 0.33
        assert call.depth_ratio == call.core_depth_ratio


@pytest.mark.bio
@pytest.mark.benchmark
def test_benchmark_df5120_real_pde_junctions():
    """Verify the 4 real benchmark PDE junctions on actual DF5120 assembly graph & ledger."""
    gfa_candidates = [
        BENCHMARK_DIR / "asm.p_utg.gfa.gz",
        BENCHMARK_DIR / "asm.p_utg.gfa",
    ]
    gfa_path = next((p for p in gfa_candidates if p.exists()), None)
    ledger_path = BENCHMARK_DIR / "df5120_end_evidence_ledger.tsv"

    if gfa_path is None or not ledger_path.exists():
        pytest.skip(f"Benchmark DF5120 dataset not present in {BENCHMARK_DIR}")

    graph = load_gfa(gfa_path)
    ledger = _load_benchmark_ledger(ledger_path)

    expected_junctions = [
        # Megabase chromosome arm ends
        {
            "port_str": "utg000001l+",
            "name": "utg000001l:3p",
            "expected_somatic": "utg000142l-",
            "expected_germline": "utg000101l+",
        },
        {
            "port_str": "utg000001l-",
            "name": "utg000001l:5p",
            "expected_somatic": "utg000076l+",
            "expected_germline": "utg000039l+",
        },
        {
            "port_str": "utg000002l+",
            "name": "utg000002l:3p",
            "expected_somatic": "utg000194l+",
            "expected_germline": "utg000115l+",
        },
        {
            "port_str": "utg000002l-",
            "name": "utg000002l:5p",
            "expected_somatic": "utg000093l-",
            "expected_germline": "utg000103l+",
        },
        {
            "port_str": "utg000009l-",
            "name": "utg000009l:5p",
            "expected_somatic": "utg000211l-",
            "expected_germline": "utg000069l+",
        },
        {
            "port_str": "utg000013l+",
            "name": "utg000013l:3p",
            "expected_somatic": "utg000167l+",
            "expected_germline": "utg000033l-",
        },
        {
            "port_str": "utg000027l-",
            "name": "utg000027l:5p",
            "expected_somatic": "utg000154l-",
            "expected_germline": "utg000189l-",
        },
        # Chr X internal breaks
        {
            "port_str": "utg000016l-",
            "name": "utg000016l:5p",
            "expected_somatic": "utg000088l+",
            "expected_germline": "utg000096l+",
        },
        {
            "port_str": "utg000019l+",
            "name": "utg000019l:3p",
            "expected_somatic": "utg000200l-",
            "expected_germline": "utg000071l+",
        },
        {
            "port_str": "utg000019l-",
            "name": "utg000019l:5p",
            "expected_somatic": "utg000129l-",
            "expected_germline": "utg000145l+",
        },
        {
            "port_str": "utg000021l-",
            "name": "utg000021l:5p",
            "expected_somatic": "utg000051l+",
            "expected_germline": "utg000071l-",
        },
        {
            "port_str": "utg000027l+",
            "name": "utg000027l:3p",
            "expected_somatic": "utg000209l+",
            "expected_germline": "utg000145l-",
        },
    ]

    for item in expected_junctions:
        port = graph._resolve_port(item["port_str"])
        call: PdeJunctionCall = classify_pde_junction(port, graph, ledger)

        assert call.is_pde_junction is True, f"Failed to identify PDE junction {item['name']}"
        assert call.port == port
        assert call.somatic_port is not None
        assert call.somatic_port.to_symbol() == item["expected_somatic"]
        assert call.germline_port is not None
        assert 0.04 <= call.depth_ratio <= 0.33
        assert 0.04 <= call.core_depth_ratio <= 0.33
        assert call.depth_ratio == call.core_depth_ratio

    # Negative controls: rDNA fork and assembly bubble must NOT be called as PDE junctions
    negative_controls = [
        ("utg000020l+", "rDNA fork"),
        ("utg000065l-", "assembly bubble"),
    ]
    for neg_port_str, neg_desc in negative_controls:
        neg_port = graph._resolve_port(neg_port_str)
        neg_call: PdeJunctionCall = classify_pde_junction(neg_port, graph, ledger)
        assert neg_call.is_pde_junction is False, (
            f"Negative control {neg_port_str} ({neg_desc}) was falsely called as PDE junction"
        )
