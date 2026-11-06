"""Unit tests for 3-axis port classifier and 1-hop PDE junction motif detector."""

from gfapdefinisher.classifier import (
    PdeJunctionCall,
    PortClassification,
    Provenance,
    TerminalFeature,
    Topology,
    classify_pde_junction,
    classify_port,
    compute_core_coverage,
)
from gfapdefinisher.graph import BidirectedGraph
from gfapdefinisher.ledger import EndEvidence
from gfapdefinisher.models import SegmentPort, Strand
from gfapdefinisher.policy import Policy, PolicyValue, StoichiometryPolicy, TelomerePolicy


def _make_dummy_evidence(
    segment_name: str,
    end: str,
    depth: float,
    degree: int = 0,
    is_tip: bool = True,
    has_terminal_telomere: bool = False,
    has_terminal_rdna: bool = False,
    usco_hit_count: int = 0,
    tiara_class: str = "",
) -> EndEvidence:
    strand = Strand.MINUS if end == "5p" else Strand.PLUS
    port = SegmentPort(segment_name, strand)
    return EndEvidence(
        segment_name=segment_name,
        end=end,
        port=port,
        length=10000,
        depth=depth,
        degree=degree,
        is_tip=is_tip,
        outgoing_neighbors=(),
        has_terminal_telomere=has_terminal_telomere,
        has_terminal_rdna=has_terminal_rdna,
        usco_hit_count=usco_hit_count,
        tiara_class=tiara_class,
    )


# ==============================================================================
# 3-Axis Port Classifier Tests
# ==============================================================================


def test_canonical_t2t_tip_classification():
    """Verify canonical T2T tip port classifies as nuclear, tip, telomeric_cap."""
    graph = BidirectedGraph()
    graph.add_segment("chr1_tip", length=50000, depth=60.0)
    port = SegmentPort("chr1_tip", Strand.FORWARD)

    ev = _make_dummy_evidence(
        "chr1_tip",
        end="3p",
        depth=60.0,
        degree=0,
        is_tip=True,
        has_terminal_telomere=True,
        usco_hit_count=5,
        tiara_class="eukarya",
    )

    call: PortClassification = classify_port(port, graph, ev)
    assert call.port == port
    assert call.provenance == Provenance.NUCLEAR
    assert call.topology == Topology.TIP
    assert call.terminal_feature == TerminalFeature.TELOMERIC_CAP


def test_collinear_internal_extension_classification():
    """Verify internal collinear port classifies as nuclear, collinear, none."""
    graph = BidirectedGraph()
    graph.add_segment("segA", length=100000, depth=60.0)
    graph.add_segment("segB", length=100000, depth=60.0)
    pA = SegmentPort("segA", Strand.FORWARD)
    pB = SegmentPort("segB", Strand.FORWARD)
    graph.add_edge(pA, pB, cigar="1000M")

    evA = _make_dummy_evidence(
        "segA",
        end="3p",
        depth=60.0,
        degree=1,
        is_tip=False,
        has_terminal_telomere=False,
        usco_hit_count=12,
        tiara_class="eukarya",
    )

    call = classify_port(pA, graph, evA)
    assert call.provenance == Provenance.NUCLEAR
    assert call.topology == Topology.COLLINEAR
    assert call.terminal_feature == TerminalFeature.NONE


def test_rdna_tangle_classification():
    """Verify multivalent rDNA array port classifies as nuclear, fork/complex, rdna_cluster."""
    graph = BidirectedGraph()
    graph.add_segment("rdna_entry", length=20000, depth=120.0)
    graph.add_segment("rdna_rep1", length=10000, depth=120.0)
    graph.add_segment("rdna_rep2", length=10000, depth=120.0)
    p_entry = SegmentPort("rdna_entry", Strand.FORWARD)
    graph.add_edge(p_entry, SegmentPort("rdna_rep1", Strand.FORWARD), cigar="500M")
    graph.add_edge(p_entry, SegmentPort("rdna_rep2", Strand.FORWARD), cigar="500M")

    ev = _make_dummy_evidence(
        "rdna_entry",
        end="3p",
        depth=120.0,
        degree=2,
        is_tip=False,
        has_terminal_rdna=True,
        usco_hit_count=0,
        tiara_class="eukarya",
    )

    call = classify_port(p_entry, graph, ev)
    assert call.provenance == Provenance.NUCLEAR
    assert call.topology == Topology.FORK
    assert call.terminal_feature == TerminalFeature.RDNA_CLUSTER


def test_organelle_and_contaminant_provenance():
    """Verify Tiara organelle and bacterial calls map to respective provenances."""
    graph = BidirectedGraph()
    graph.add_segment("mito1", length=15000, depth=300.0)
    graph.add_segment("bact1", length=50000, depth=20.0)

    ev_mito = _make_dummy_evidence("mito1", "3p", depth=300.0, tiara_class="organelle")
    call_mito = classify_port(SegmentPort("mito1", Strand.FORWARD), graph, ev_mito)
    assert call_mito.provenance == Provenance.ORGANELLE

    ev_bact = _make_dummy_evidence("bact1", "5p", depth=20.0, tiara_class="bacteria")
    call_bact = classify_port(SegmentPort("bact1", Strand.MINUS), graph, ev_bact)
    assert call_bact.provenance == Provenance.CONTAMINANT

    ev_unk = _make_dummy_evidence("unk1", "3p", depth=5.0, tiara_class="unknown")
    call_unk = classify_port(SegmentPort("unk1", Strand.FORWARD), graph, ev_unk)
    assert call_unk.provenance == Provenance.UNKNOWN


# ==============================================================================
# PDE Junction Motif Detector Tests
# ==============================================================================


def test_canonical_pde_junction_detection_direct_cap():
    """Verify synthetic 60x stem branching to 50x telomere cap and 6x GRS branch."""
    graph = BidirectedGraph()
    graph.add_segment("stem", length=100000, depth=60.0)
    graph.add_segment("somatic_cap", length=25000, depth=50.0)
    graph.add_segment("germline_grs", length=80000, depth=6.0)

    p_stem = SegmentPort("stem", Strand.FORWARD)
    p_som = SegmentPort("somatic_cap", Strand.FORWARD)
    p_grs = SegmentPort("germline_grs", Strand.FORWARD)

    graph.add_edge(p_stem, p_som, cigar="100M")
    graph.add_edge(p_stem, p_grs, cigar="100M")

    ev_stem = _make_dummy_evidence("stem", "3p", depth=60.0, degree=2, is_tip=False, usco_hit_count=10)
    ev_som = _make_dummy_evidence("somatic_cap", "3p", depth=50.0, degree=0, is_tip=True, has_terminal_telomere=True)
    ev_grs = _make_dummy_evidence("germline_grs", "3p", depth=6.0, degree=0, is_tip=True, has_terminal_telomere=False)

    ledger = {
        p_stem.to_symbol(): ev_stem,
        p_som.to_symbol(): ev_som,
        p_grs.to_symbol(): ev_grs,
    }

    call: PdeJunctionCall = classify_pde_junction(p_stem, graph, ledger)
    assert call.is_pde_junction is True
    assert call.port == p_stem
    assert call.somatic_port == p_som
    assert call.germline_port == p_grs
    assert abs(call.depth_ratio - 0.10) < 1e-4


def test_1hop_pde_junction_detection_linked_cap():
    """Verify 1-hop reachability: stem -> intermediate unitig (55x) -> telomere cap (50x)."""
    graph = BidirectedGraph()
    graph.add_segment("stem", length=100000, depth=60.0)
    graph.add_segment("som_mid", length=20000, depth=55.0)
    graph.add_segment("som_cap", length=15000, depth=50.0)
    graph.add_segment("germline_grs", length=80000, depth=6.0)

    p_stem = SegmentPort("stem", Strand.FORWARD)
    p_som_mid = SegmentPort("som_mid", Strand.FORWARD)
    p_som_cap = SegmentPort("som_cap", Strand.FORWARD)
    p_grs = SegmentPort("germline_grs", Strand.FORWARD)

    # Stem bifurcates into som_mid and germline_grs
    graph.add_edge(p_stem, p_som_mid, cigar="100M")
    graph.add_edge(p_stem, p_grs, cigar="100M")
    # som_mid links 1 hop further to som_cap (telomere cap tip)
    graph.add_edge(p_som_mid, p_som_cap, cigar="100M")

    ev_stem = _make_dummy_evidence("stem", "3p", depth=60.0, degree=2, is_tip=False)
    ev_som_mid = _make_dummy_evidence("som_mid", "3p", depth=55.0, degree=1, is_tip=False)
    ev_som_cap = _make_dummy_evidence("som_cap", "3p", depth=50.0, degree=0, is_tip=True, has_terminal_telomere=True)
    ev_grs = _make_dummy_evidence("germline_grs", "3p", depth=6.0, degree=0, is_tip=True)

    ledger = {
        p_stem.to_symbol(): ev_stem,
        p_som_mid.to_symbol(): ev_som_mid,
        p_som_cap.to_symbol(): ev_som_cap,
        p_grs.to_symbol(): ev_grs,
    }

    call = classify_pde_junction(p_stem, graph, ledger)
    assert call.is_pde_junction is True
    assert call.somatic_port == p_som_mid
    assert call.germline_port == p_grs
    assert abs(call.depth_ratio - 0.10) < 1e-4


def test_negative_cases_symmetric_bubble():
    """Verify heterozygous bubble (two 30x non-telomeric branches) is rejected."""
    graph = BidirectedGraph()
    graph.add_segment("stem", length=100000, depth=60.0)
    graph.add_segment("branch1", length=50000, depth=30.0)
    graph.add_segment("branch2", length=50000, depth=30.0)

    p_stem = SegmentPort("stem", Strand.FORWARD)
    p_b1 = SegmentPort("branch1", Strand.FORWARD)
    p_b2 = SegmentPort("branch2", Strand.FORWARD)

    graph.add_edge(p_stem, p_b1, cigar="100M")
    graph.add_edge(p_stem, p_b2, cigar="100M")

    ev_stem = _make_dummy_evidence("stem", "3p", depth=60.0, degree=2, is_tip=False)
    ev_b1 = _make_dummy_evidence("branch1", "3p", depth=30.0, degree=0, is_tip=True)
    ev_b2 = _make_dummy_evidence("branch2", "3p", depth=30.0, degree=0, is_tip=True)

    ledger = {
        p_stem.to_symbol(): ev_stem,
        p_b1.to_symbol(): ev_b1,
        p_b2.to_symbol(): ev_b2,
    }

    call = classify_pde_junction(p_stem, graph, ledger)
    assert call.is_pde_junction is False


def test_negative_cases_depth_ratio_too_low():
    """Verify germline depth ratio below 0.04 (e.g. 1x / 60x = 0.0167) is rejected."""
    graph = BidirectedGraph()
    graph.add_segment("stem", length=100000, depth=60.0)
    graph.add_segment("somatic_cap", length=25000, depth=50.0)
    graph.add_segment("noise_branch", length=5000, depth=1.0)

    p_stem = SegmentPort("stem", Strand.FORWARD)
    p_som = SegmentPort("somatic_cap", Strand.FORWARD)
    p_noise = SegmentPort("noise_branch", Strand.FORWARD)

    graph.add_edge(p_stem, p_som, cigar="100M")
    graph.add_edge(p_stem, p_noise, cigar="100M")

    ledger = {
        p_stem.to_symbol(): _make_dummy_evidence("stem", "3p", depth=60.0, degree=2),
        p_som.to_symbol(): _make_dummy_evidence("somatic_cap", "3p", depth=50.0, has_terminal_telomere=True),
        p_noise.to_symbol(): _make_dummy_evidence("noise_branch", "3p", depth=1.0),
    }

    call = classify_pde_junction(p_stem, graph, ledger)
    assert call.is_pde_junction is False


def test_negative_cases_depth_ratio_too_high():
    """Verify germline candidate ratio above 0.40 (e.g. 30x / 60x = 0.50) is rejected."""
    graph = BidirectedGraph()
    graph.add_segment("stem", length=100000, depth=60.0)
    graph.add_segment("somatic_cap", length=25000, depth=50.0)
    graph.add_segment("high_branch", length=50000, depth=30.0)

    p_stem = SegmentPort("stem", Strand.FORWARD)
    p_som = SegmentPort("somatic_cap", Strand.FORWARD)
    p_high = SegmentPort("high_branch", Strand.FORWARD)

    graph.add_edge(p_stem, p_som, cigar="100M")
    graph.add_edge(p_stem, p_high, cigar="100M")

    ledger = {
        p_stem.to_symbol(): _make_dummy_evidence("stem", "3p", depth=60.0, degree=2),
        p_som.to_symbol(): _make_dummy_evidence("somatic_cap", "3p", depth=50.0, has_terminal_telomere=True),
        p_high.to_symbol(): _make_dummy_evidence("high_branch", "3p", depth=30.0),
    }

    call = classify_pde_junction(p_stem, graph, ledger)
    assert call.is_pde_junction is False


def test_negative_cases_missing_telomere_cap():
    """Verify somatic candidate lacking telomeric cap is rejected."""
    graph = BidirectedGraph()
    graph.add_segment("stem", length=100000, depth=60.0)
    graph.add_segment("som_uncapped", length=25000, depth=50.0)
    graph.add_segment("germline_grs", length=80000, depth=6.0)

    p_stem = SegmentPort("stem", Strand.FORWARD)
    p_som = SegmentPort("som_uncapped", Strand.FORWARD)
    p_grs = SegmentPort("germline_grs", Strand.FORWARD)

    graph.add_edge(p_stem, p_som, cigar="100M")
    graph.add_edge(p_stem, p_grs, cigar="100M")

    ledger = {
        p_stem.to_symbol(): _make_dummy_evidence("stem", "3p", depth=60.0, degree=2),
        p_som.to_symbol(): _make_dummy_evidence("som_uncapped", "3p", depth=50.0, has_terminal_telomere=False),
        p_grs.to_symbol(): _make_dummy_evidence("germline_grs", "3p", depth=6.0),
    }

    call = classify_pde_junction(p_stem, graph, ledger)
    assert call.is_pde_junction is False


def test_negative_cases_linear_tip_or_zero_depth():
    """Verify degree < 2 and zero depth stems are rejected."""
    graph = BidirectedGraph()
    graph.add_segment("stem", length=100000, depth=0.0)
    p_stem = SegmentPort("stem", Strand.FORWARD)

    call_zero = classify_pde_junction(p_stem, graph, {})
    assert call_zero.is_pde_junction is False


def test_custom_policy_override():
    """Verify custom stoichiometry policy thresholds are respected."""
    graph = BidirectedGraph()
    graph.add_segment("stem", length=100000, depth=100.0)
    graph.add_segment("som_cap", length=25000, depth=90.0)
    graph.add_segment("grs", length=50000, depth=45.0)

    p_stem = SegmentPort("stem", Strand.FORWARD)
    p_som = SegmentPort("som_cap", Strand.FORWARD)
    p_grs = SegmentPort("grs", Strand.FORWARD)

    graph.add_edge(p_stem, p_som, cigar="100M")
    graph.add_edge(p_stem, p_grs, cigar="100M")

    ledger = {
        p_stem.to_symbol(): _make_dummy_evidence("stem", "3p", depth=100.0, degree=2),
        p_som.to_symbol(): _make_dummy_evidence("som_cap", "3p", depth=90.0, has_terminal_telomere=True),
        p_grs.to_symbol(): _make_dummy_evidence("grs", "3p", depth=45.0),
    }

    # Under default policy (max_germline_depth_ratio=0.40), 45/100 = 0.45 is too high
    call_default = classify_pde_junction(p_stem, graph, ledger)
    assert call_default.is_pde_junction is False

    # Under custom policy allowing max_ratio=0.50, it passes
    custom_policy = Policy(
        telomere=TelomerePolicy(),
        stoichiometry=StoichiometryPolicy(
            min_germline_depth_ratio=PolicyValue(value=0.05, justification="Custom test", provenance="Test"),
            max_germline_depth_ratio=PolicyValue(value=0.50, justification="Custom test", provenance="Test"),
        ),
    )
    call_custom = classify_pde_junction(p_stem, graph, ledger, policy=custom_policy)
    assert call_custom.is_pde_junction is True
    assert call_custom.somatic_port == p_som
    assert call_custom.germline_port == p_grs
    assert abs(call_custom.depth_ratio - 0.45) < 1e-4


# ==============================================================================
# Core Coverage, Provenance Extensions, and Structural Guardrail Unit Tests
# ==============================================================================


def test_compute_core_coverage_tier1_usco_and_eukarya():
    """Verify compute_core_coverage calculates length-weighted median for Tier 1 segments."""
    graph = BidirectedGraph()
    graph.add_segment("u1", length=1000, depth=20.0)
    graph.add_segment("u2", length=3000, depth=60.0)
    graph.add_segment("u3", length=5000, depth=100.0)

    p1 = SegmentPort("u1", Strand.FORWARD)
    p2 = SegmentPort("u2", Strand.FORWARD)
    p3 = SegmentPort("u3", Strand.FORWARD)

    ledger = {
        p1.to_symbol(): _make_dummy_evidence("u1", "3p", depth=20.0, usco_hit_count=3, tiara_class="eukarya"),
        p2.to_symbol(): _make_dummy_evidence("u2", "3p", depth=60.0, usco_hit_count=10, tiara_class="eukarya"),
        p3.to_symbol(): _make_dummy_evidence("u3", "3p", depth=100.0, usco_hit_count=0, tiara_class="bacteria"),
    }

    # Eligible pool: u1 (1000 bp @ 20x) and u2 (3000 bp @ 60x). Total: 4000 bp.
    # Half-point: 2000 bp -> falls into u2 (cum_len 4000 >= 2000). Length-weighted median = 60.0.
    core_cov = compute_core_coverage(graph, ledger)
    assert abs(core_cov - 60.0) < 1e-4


def test_compute_core_coverage_fallback_tiers():
    """Verify compute_core_coverage falls back to Tier 2 (eukarya) and Tier 3 (all segments)."""
    # Tier 2 fallback: no USCO hits, but eukarya present
    g2 = BidirectedGraph()
    g2.add_segment("s1", length=1000, depth=30.0)
    g2.add_segment("s2", length=2000, depth=50.0)
    p1 = SegmentPort("s1", Strand.FORWARD)
    p2 = SegmentPort("s2", Strand.FORWARD)
    l2 = {
        p1.to_symbol(): _make_dummy_evidence("s1", "3p", depth=30.0, usco_hit_count=0, tiara_class="eukarya"),
        p2.to_symbol(): _make_dummy_evidence("s2", "3p", depth=50.0, usco_hit_count=0, tiara_class="eukarya"),
    }
    # Total: 3000 bp. Half: 1500 bp -> s2 (cum 3000 >= 1500) -> 50.0
    assert abs(compute_core_coverage(g2, l2) - 50.0) < 1e-4

    # Tier 3 fallback: no USCO, no eukarya -> all segments
    g3 = BidirectedGraph()
    g3.add_segment("b1", length=1000, depth=15.0)
    g3.add_segment("b2", length=2000, depth=40.0)
    p_b1 = SegmentPort("b1", Strand.FORWARD)
    p_b2 = SegmentPort("b2", Strand.FORWARD)
    l3 = {
        p_b1.to_symbol(): _make_dummy_evidence("b1", "3p", depth=15.0, tiara_class="unknown"),
        p_b2.to_symbol(): _make_dummy_evidence("b2", "3p", depth=40.0, tiara_class="unknown"),
    }
    assert abs(compute_core_coverage(g3, l3) - 40.0) < 1e-4

    # Empty graph
    g_empty = BidirectedGraph()
    assert compute_core_coverage(g_empty, {}) == 0.0


def test_prokarya_and_mitochondria_provenance():
    """Verify mitochondria and prokarya Tiara classes map to ORGANELLE and CONTAMINANT."""
    graph = BidirectedGraph()
    graph.add_segment("mito_seg", length=14000, depth=250.0)
    graph.add_segment("prok_seg", length=40000, depth=15.0)

    ev_mito = _make_dummy_evidence("mito_seg", "3p", depth=250.0, tiara_class="mitochondria")
    call_mito = classify_port(SegmentPort("mito_seg", Strand.FORWARD), graph, ev_mito)
    assert call_mito.provenance == Provenance.ORGANELLE

    ev_prok = _make_dummy_evidence("prok_seg", "5p", depth=15.0, tiara_class="prokarya")
    call_prok = classify_port(SegmentPort("prok_seg", Strand.MINUS), graph, ev_prok)
    assert call_prok.provenance == Provenance.CONTAMINANT


def test_stem_rdna_rejection():
    """Verify stem carrying rDNA is rejected as a PDE junction (prevents rDNA fork false positives)."""
    graph = BidirectedGraph()
    graph.add_segment("rdna_stem", length=40000, depth=60.0)
    graph.add_segment("som_cap", length=20000, depth=50.0)
    graph.add_segment("grs_branch", length=50000, depth=6.0)

    p_stem = SegmentPort("rdna_stem", Strand.FORWARD)
    p_som = SegmentPort("som_cap", Strand.FORWARD)
    p_grs = SegmentPort("grs_branch", Strand.FORWARD)

    graph.add_edge(p_stem, p_som, cigar="100M")
    graph.add_edge(p_stem, p_grs, cigar="100M")

    # Stem with rdna_hit_count > 0
    ev_stem_hits = _make_dummy_evidence(
        "rdna_stem", "3p", depth=60.0, degree=2, is_tip=False, has_terminal_rdna=False, usco_hit_count=0
    )
    ev_stem_hits = EndEvidence(
        segment_name="rdna_stem",
        end="3p",
        port=p_stem,
        length=40000,
        depth=60.0,
        degree=2,
        is_tip=False,
        outgoing_neighbors=(),
        has_terminal_telomere=False,
        has_terminal_rdna=False,
        rdna_hit_count=4,
        usco_hit_count=0,
        tiara_class="eukarya",
    )
    ev_som = _make_dummy_evidence("som_cap", "3p", depth=50.0, degree=0, is_tip=True, has_terminal_telomere=True)
    ev_grs = _make_dummy_evidence("grs_branch", "3p", depth=6.0, degree=0, is_tip=True)

    ledger_hits = {
        p_stem.to_symbol(): ev_stem_hits,
        p_som.to_symbol(): ev_som,
        p_grs.to_symbol(): ev_grs,
    }
    call_hits = classify_pde_junction(p_stem, graph, ledger_hits)
    assert call_hits.is_pde_junction is False

    # Stem with has_terminal_rdna = True
    ev_stem_term = _make_dummy_evidence("rdna_stem", "3p", depth=60.0, degree=2, is_tip=False, has_terminal_rdna=True)
    ledger_term = {
        p_stem.to_symbol(): ev_stem_term,
        p_som.to_symbol(): ev_som,
        p_grs.to_symbol(): ev_grs,
    }
    call_term = classify_pde_junction(p_stem, graph, ledger_term)
    assert call_term.is_pde_junction is False


def test_rejoining_bubble_rejection():
    """Verify bifurcations reconverging onto the same downstream target segment are rejected."""
    graph = BidirectedGraph()
    graph.add_segment("bubble_stem", length=70000, depth=55.0)
    graph.add_segment("arm_a", length=30000, depth=28.0)
    graph.add_segment("arm_b", length=25000, depth=18.0)
    graph.add_segment("target_seg", length=80000, depth=55.0)

    p_stem = SegmentPort("bubble_stem", Strand.FORWARD)
    p_arm_a = SegmentPort("arm_a", Strand.FORWARD)
    p_arm_b = SegmentPort("arm_b", Strand.FORWARD)
    p_target = SegmentPort("target_seg", Strand.FORWARD)

    # Stem forks into arm_a and arm_b
    graph.add_edge(p_stem, p_arm_a, cigar="100M")
    graph.add_edge(p_stem, p_arm_b, cigar="100M")
    # Both arms reconverge onto target_seg
    graph.add_edge(p_arm_a, p_target, cigar="100M")
    graph.add_edge(p_arm_b, p_target, cigar="100M")

    ledger = {
        p_stem.to_symbol(): _make_dummy_evidence("bubble_stem", "3p", depth=55.0, degree=2, is_tip=False),
        p_arm_a.to_symbol(): _make_dummy_evidence("arm_a", "3p", depth=28.0, degree=1, is_tip=False),
        p_arm_b.to_symbol(): _make_dummy_evidence("arm_b", "3p", depth=18.0, degree=1, is_tip=False),
        p_target.to_symbol(): _make_dummy_evidence("target_seg", "5p", depth=55.0, degree=0, is_tip=True),
    }

    call = classify_pde_junction(p_stem, graph, ledger)
    assert call.is_pde_junction is False


def test_complex_grs_repeat_entry():
    """Verify complex GRS repeat entry is identified when somatic cap is verified but no single-copy GRS."""
    graph = BidirectedGraph()
    graph.add_segment("stem", length=100000, depth=60.0)
    graph.add_segment("somatic_cap", length=25000, depth=50.0)
    graph.add_segment("repeat_grs", length=15000, depth=60.0)  # Multi-copy depth, degree 2
    graph.add_segment("satellite1", length=5000, depth=60.0)
    graph.add_segment("satellite2", length=5000, depth=60.0)

    p_stem = SegmentPort("stem", Strand.FORWARD)
    p_som = SegmentPort("somatic_cap", Strand.FORWARD)
    p_rep = SegmentPort("repeat_grs", Strand.FORWARD)

    graph.add_edge(p_stem, p_som, cigar="100M")
    graph.add_edge(p_stem, p_rep, cigar="100M")
    # repeat_grs branches into satellites (degree >= 2), does not loop back to stem
    graph.add_edge(p_rep, SegmentPort("satellite1", Strand.FORWARD), cigar="50M")
    graph.add_edge(p_rep, SegmentPort("satellite2", Strand.FORWARD), cigar="50M")

    ledger = {
        p_stem.to_symbol(): _make_dummy_evidence("stem", "3p", depth=60.0, degree=2, is_tip=False, usco_hit_count=5),
        p_som.to_symbol(): _make_dummy_evidence(
            "somatic_cap", "3p", depth=50.0, degree=0, is_tip=True, has_terminal_telomere=True
        ),
        p_rep.to_symbol(): _make_dummy_evidence("repeat_grs", "3p", depth=60.0, degree=2, is_tip=False),
    }

    call = classify_pde_junction(p_stem, graph, ledger)
    assert call.is_pde_junction is True
    assert call.somatic_port == p_som
    assert call.germline_port == p_rep
    assert call.is_repeat_entry is True
