"""Unit tests for 3-axis port classifier and 1-hop PDE junction motif detector."""

from gfapdefinisher.classifier import (
    PdeJunctionCall,
    PortClassification,
    Provenance,
    TerminalFeature,
    Topology,
    classify_pde_junction,
    classify_port,
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
    assert call.confidence == "HIGH"


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
