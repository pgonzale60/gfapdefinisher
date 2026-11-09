"""Synthetic biological falsification test suite for BIO-001 through BIO-006.

Validates that assembly graph algorithms adhere strictly to biological priors
defined in docs/biology/RULES.md:
- BIO-001: Unconditional Terminal Completion (rDNA rejection; telomere requirement)
- BIO-002: Composition & Stoichiometry Over Arbitrary Length Cutoffs (micro-GRS <5 kb preserved)
- BIO-003: Graph-Inferred PDE Breakpoints (coordinate-free topological discovery)
- BIO-004: Repeat Representation Preservation for Variable Satellites (representative sequence retained)
- BIO-005: Lossless Overlap and Provenance Invariance (exact CIGAR, tags, and A-lines)
- BIO-006: Held-Out Generalization Guard (permutation and naming order invariance)
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

from gfapdefinisher.cigar import cigar_overlap_length
from gfapdefinisher.classifier import (
    PdeJunctionCall,
    PortClassification,
    TerminalFeature,
    classify_pde_junction,
    classify_port,
)
from gfapdefinisher.graph import BidirectedGraph, load_gfa
from gfapdefinisher.ledger import EndEvidence
from gfapdefinisher.models import SegmentPort, Strand
from gfapdefinisher.policy import DEFAULT_POLICY


def _make_evidence(
    segment_name: str,
    end: str,
    depth: float,
    length: int = 50000,
    degree: int = 1,
    is_tip: bool = False,
    has_terminal_telomere: bool = False,
    has_terminal_rdna: bool = False,
    rdna_hit_count: int = 0,
    usco_hit_count: int = 0,
    tiara_class: str = "eukarya",
    outgoing_neighbors: tuple[str, ...] = (),
) -> EndEvidence:
    """Build synthetic EndEvidence fixture for tests."""
    strand = Strand.MINUS if end == "5p" else Strand.PLUS
    port = SegmentPort(segment_name, strand)
    return EndEvidence(
        segment_name=segment_name,
        end=end,
        port=port,
        length=length,
        depth=depth,
        degree=degree,
        is_tip=is_tip,
        outgoing_neighbors=outgoing_neighbors,
        telomere_canonical_count=50 if has_terminal_telomere else 0,
        telomere_opposite_count=0,
        has_terminal_telomere=has_terminal_telomere,
        rdna_hit_count=rdna_hit_count,
        has_terminal_rdna=has_terminal_rdna,
        usco_hit_count=usco_hit_count,
        tiara_class=tiara_class,
    )


# ==============================================================================
# BIO-001: Unconditional Terminal Completion
# ==============================================================================


@pytest.mark.bio
def test_bio_001_rdna_array_rejection_requires_telomeric_cap():
    """BIO-001: Chromosome arms entering rDNA arrays must not be treated as legitimate termini.

    An authentic chromosome terminus requires canonical telomeric repeats; naked rDNA
    or satellite clusters are internal repeats, not telomeric caps.
    """
    graph = BidirectedGraph()
    # Chromosome arm entering an rDNA cassette
    graph.add_segment("chr_arm", length=120_000, depth=62.0)
    graph.add_segment("rdna_cassette", length=18_000, depth=240.0)
    graph.add_segment("grs_branch", length=40_000, depth=7.0)
    graph.add_segment("tel_somatic_cap", length=25_000, depth=45.0)

    p_arm = SegmentPort("chr_arm", Strand.FORWARD)
    p_rdna = SegmentPort("rdna_cassette", Strand.FORWARD)
    p_grs = SegmentPort("grs_branch", Strand.FORWARD)
    p_tel = SegmentPort("tel_somatic_cap", Strand.FORWARD)

    # Fork 1: Arm enters rDNA cluster and GRS
    graph.add_edge(p_arm, p_rdna, cigar="100M")
    graph.add_edge(p_arm, p_grs, cigar="100M")

    ev_arm = _make_evidence("chr_arm", "3p", depth=62.0, degree=2, is_tip=False, usco_hit_count=10)
    ev_rdna = _make_evidence(
        "rdna_cassette",
        "3p",
        depth=240.0,
        degree=0,
        is_tip=True,
        has_terminal_rdna=True,
        rdna_hit_count=15,
        has_terminal_telomere=False,
    )
    ev_grs = _make_evidence("grs_branch", "3p", depth=7.0, degree=0, is_tip=True)

    ledger = {
        p_arm.to_symbol(): ev_arm,
        p_rdna.to_symbol(): ev_rdna,
        p_grs.to_symbol(): ev_grs,
    }

    # 1. Port classification confirms rDNA is recognized as RDNA_CLUSTER, not TELOMERIC_CAP
    classif: PortClassification = classify_port(p_rdna, graph, ev_rdna)
    assert classif.terminal_feature == TerminalFeature.RDNA_CLUSTER
    assert classif.terminal_feature != TerminalFeature.TELOMERIC_CAP

    # 2. PDE junction motif MUST reject fork entering rDNA as a somatic terminus
    call_rdna: PdeJunctionCall = classify_pde_junction(p_arm, graph, ledger)
    assert call_rdna.is_pde_junction is False, "Fork entering rDNA must not be called as a somatic diminution cap"

    # 3. Now connect the arm to an authentic telomeric cap
    graph_tel = BidirectedGraph()
    graph_tel.add_segment("chr_arm", length=120_000, depth=62.0)
    graph_tel.add_segment("tel_somatic_cap", length=25_000, depth=45.0)
    graph_tel.add_segment("grs_branch", length=40_000, depth=7.0)
    graph_tel.add_edge(p_arm, p_tel, cigar="100M")
    graph_tel.add_edge(p_arm, p_grs, cigar="100M")

    ev_tel = _make_evidence("tel_somatic_cap", "3p", depth=45.0, degree=0, is_tip=True, has_terminal_telomere=True)
    ledger_tel = {
        p_arm.to_symbol(): ev_arm,
        p_tel.to_symbol(): ev_tel,
        p_grs.to_symbol(): ev_grs,
    }

    call_tel: PdeJunctionCall = classify_pde_junction(p_arm, graph_tel, ledger_tel)
    assert call_tel.is_pde_junction is True
    assert call_tel.somatic_port == p_tel
    assert call_tel.germline_port == p_grs


# ==============================================================================
# BIO-002: Composition & Stoichiometry Over Arbitrary Length Cutoffs
# ==============================================================================


@pytest.mark.bio
def test_bio_002_micro_grs_unitigs_under_5kb_not_filtered_by_length():
    """BIO-002: Micro-GRS/short unitigs (<5 kb) are not filtered out when stoichiometric conditions are satisfied.

    Decisions must be driven by physical sequence composition and quantitative depth
    stoichiometry, never by arbitrary sequence length cutoffs.
    """
    graph = BidirectedGraph()
    graph.add_segment("stem", length=100_000, depth=64.0)
    graph.add_segment("som_cap", length=25_000, depth=48.0)
    # Micro-GRS unitig with length 1,200 bp (< 5,000 bp)
    graph.add_segment("micro_grs_1200bp", length=1200, depth=7.0)

    p_stem = SegmentPort("stem", Strand.FORWARD)
    p_som = SegmentPort("som_cap", Strand.FORWARD)
    p_grs = SegmentPort("micro_grs_1200bp", Strand.FORWARD)

    graph.add_edge(p_stem, p_som, cigar="100M")
    graph.add_edge(p_stem, p_grs, cigar="100M")

    ev_stem = _make_evidence("stem", "3p", depth=64.0, degree=2, is_tip=False, usco_hit_count=5)
    ev_som = _make_evidence("som_cap", "3p", depth=48.0, degree=0, is_tip=True, has_terminal_telomere=True)
    ev_grs = _make_evidence("micro_grs_1200bp", "3p", depth=7.0, length=1200, degree=0, is_tip=True)

    ledger = {
        p_stem.to_symbol(): ev_stem,
        p_som.to_symbol(): ev_som,
        p_grs.to_symbol(): ev_grs,
    }

    call: PdeJunctionCall = classify_pde_junction(p_stem, graph, ledger)
    assert call.is_pde_junction is True
    assert call.somatic_port == p_som
    assert call.germline_port == p_grs
    assert abs(call.depth_ratio - (7.0 / 64.0)) < 1e-4

    # Also test an ultra-short unitig: 450 bp (< 1 kb)
    graph_sub_kb = BidirectedGraph()
    graph_sub_kb.add_segment("stem", length=100_000, depth=64.0)
    graph_sub_kb.add_segment("som_cap", length=25_000, depth=48.0)
    graph_sub_kb.add_segment("micro_grs_450bp", length=450, depth=6.0)

    p_grs_450 = SegmentPort("micro_grs_450bp", Strand.FORWARD)
    graph_sub_kb.add_edge(p_stem, p_som, cigar="100M")
    graph_sub_kb.add_edge(p_stem, p_grs_450, cigar="100M")

    ev_grs_450 = _make_evidence("micro_grs_450bp", "3p", depth=6.0, length=450, degree=0, is_tip=True)
    ledger_sub_kb = {
        p_stem.to_symbol(): ev_stem,
        p_som.to_symbol(): ev_som,
        p_grs_450.to_symbol(): ev_grs_450,
    }

    call_sub_kb: PdeJunctionCall = classify_pde_junction(p_stem, graph_sub_kb, ledger_sub_kb)
    assert call_sub_kb.is_pde_junction is True
    assert call_sub_kb.germline_port == p_grs_450

    # Falsification check: inspect AST across src/gfapdefinisher to ensure zero static length cutoffs on GRS
    src_dir = Path(__file__).resolve().parent.parent.parent / "src" / "gfapdefinisher"
    forbidden_thresholds = {500, 1000, 2000, 5000, 10000}
    for py_file in src_dir.rglob("*.py"):
        if "policy" in py_file.parts:
            continue
        tree = ast.parse(py_file.read_text(encoding="utf-8"), filename=str(py_file))
        for node in ast.walk(tree):
            if isinstance(node, ast.Compare):
                # Check both left-hand and right-hand comparison operands
                for operand in [node.left, *node.comparators]:
                    if (
                        isinstance(operand, ast.Constant)
                        and isinstance(operand.value, int)
                        and operand.value in forbidden_thresholds
                    ):
                        pytest.fail(
                            f"{py_file.name}:{node.lineno} contains forbidden static length cutoff {operand.value}"
                        )


# ==============================================================================
# BIO-003: Graph-Inferred Programmed DNA Elimination (PDE) Breakpoints
# ==============================================================================


@pytest.mark.bio
def test_bio_003_graph_topology_discovers_diminution_without_coordinates():
    """BIO-003: Diminution forks are inferred directly from graph topology and depth stoichiometry.

    Must operate without requiring external genomic coordinate BED files or reference alignments.
    """
    graph = BidirectedGraph()
    # Diminution fork
    graph.add_segment("somatic_stem", length=150_000, depth=60.0)
    graph.add_segment("som_cap", length=20_000, depth=52.0)
    graph.add_segment("germline_grs", length=90_000, depth=6.5)

    # Internal collinear segment (negative control 1)
    graph.add_segment("collinear_1", length=80_000, depth=60.0)
    graph.add_segment("collinear_2", length=70_000, depth=60.0)

    # Reconverging heterozygous bubble (negative control 2)
    graph.add_segment("bubble_stem", length=100_000, depth=60.0)
    graph.add_segment("bubble_arm_a", length=30_000, depth=30.0)
    graph.add_segment("bubble_arm_b", length=30_000, depth=30.0)
    graph.add_segment("bubble_reconverg", length=50_000, depth=60.0)

    # Spurious low-depth noise branch (negative control 3: depth ratio < min_germline_depth_ratio)
    noise_stem_depth = 60.0
    noise_branch_depth = noise_stem_depth * DEFAULT_POLICY.stoichiometry.min_germline_depth_ratio.value * 0.5
    graph.add_segment("noise_stem", length=100_000, depth=noise_stem_depth)
    graph.add_segment("noise_som_cap", length=20_000, depth=50.0)
    graph.add_segment("noise_branch", length=40_000, depth=noise_branch_depth)

    # Wire edges
    p_stem = SegmentPort("somatic_stem", Strand.FORWARD)
    p_som = SegmentPort("som_cap", Strand.FORWARD)
    p_grs = SegmentPort("germline_grs", Strand.FORWARD)
    graph.add_edge(p_stem, p_som, cigar="100M")
    graph.add_edge(p_stem, p_grs, cigar="100M")

    p_coll_1 = SegmentPort("collinear_1", Strand.FORWARD)
    p_coll_2 = SegmentPort("collinear_2", Strand.FORWARD)
    graph.add_edge(p_coll_1, p_coll_2, cigar="100M")

    p_bub_stem = SegmentPort("bubble_stem", Strand.FORWARD)
    p_bub_a = SegmentPort("bubble_arm_a", Strand.FORWARD)
    p_bub_b = SegmentPort("bubble_arm_b", Strand.FORWARD)
    p_bub_reconv = SegmentPort("bubble_reconverg", Strand.FORWARD)
    graph.add_edge(p_bub_stem, p_bub_a, cigar="100M")
    graph.add_edge(p_bub_stem, p_bub_b, cigar="100M")
    graph.add_edge(p_bub_a, p_bub_reconv, cigar="100M")
    graph.add_edge(p_bub_b, p_bub_reconv, cigar="100M")

    p_noise_stem = SegmentPort("noise_stem", Strand.FORWARD)
    p_noise_som = SegmentPort("noise_som_cap", Strand.FORWARD)
    p_noise_br = SegmentPort("noise_branch", Strand.FORWARD)
    graph.add_edge(p_noise_stem, p_noise_som, cigar="100M")
    graph.add_edge(p_noise_stem, p_noise_br, cigar="100M")

    ev_stem = _make_evidence("somatic_stem", "3p", depth=60.0, degree=2, is_tip=False, usco_hit_count=8)
    ev_som = _make_evidence("som_cap", "3p", depth=52.0, degree=0, is_tip=True, has_terminal_telomere=True)
    ev_grs = _make_evidence("germline_grs", "3p", depth=6.5, degree=0, is_tip=True)

    ev_coll_1 = _make_evidence("collinear_1", "3p", depth=60.0, degree=1, is_tip=False, usco_hit_count=5)
    ev_coll_2 = _make_evidence("collinear_2", "3p", depth=60.0, degree=1, is_tip=False, usco_hit_count=5)

    ev_bub_stem = _make_evidence("bubble_stem", "3p", depth=60.0, degree=2, is_tip=False)
    ev_bub_a = _make_evidence("bubble_arm_a", "3p", depth=30.0, degree=1, is_tip=False)
    ev_bub_b = _make_evidence("bubble_arm_b", "3p", depth=30.0, degree=1, is_tip=False)
    ev_bub_reconv = _make_evidence("bubble_reconverg", "3p", depth=60.0, degree=1, is_tip=False)

    ev_noise_stem = _make_evidence("noise_stem", "3p", depth=noise_stem_depth, degree=2, is_tip=False)
    ev_noise_som = _make_evidence("noise_som_cap", "3p", depth=50.0, degree=0, is_tip=True, has_terminal_telomere=True)
    ev_noise_br = _make_evidence("noise_branch", "3p", depth=noise_branch_depth, degree=0, is_tip=True)

    ledger = {
        p_stem.to_symbol(): ev_stem,
        p_som.to_symbol(): ev_som,
        p_grs.to_symbol(): ev_grs,
        p_coll_1.to_symbol(): ev_coll_1,
        p_coll_2.to_symbol(): ev_coll_2,
        p_bub_stem.to_symbol(): ev_bub_stem,
        p_bub_a.to_symbol(): ev_bub_a,
        p_bub_b.to_symbol(): ev_bub_b,
        p_bub_reconv.to_symbol(): ev_bub_reconv,
        p_noise_stem.to_symbol(): ev_noise_stem,
        p_noise_som.to_symbol(): ev_noise_som,
        p_noise_br.to_symbol(): ev_noise_br,
    }

    # Authentic diminution fork discovered without external coordinates
    call: PdeJunctionCall = classify_pde_junction(p_stem, graph, ledger)
    assert call.is_pde_junction is True
    assert call.somatic_port == p_som
    assert call.germline_port == p_grs

    # Negative control 1: Collinear segment is not a PDE fork (degree = 1)
    call_coll = classify_pde_junction(p_coll_1, graph, ledger)
    assert call_coll.is_pde_junction is False

    # Negative control 2: Reconverging bubble is rejected
    call_bub = classify_pde_junction(p_bub_stem, graph, ledger)
    assert call_bub.is_pde_junction is False

    # Negative control 3: Spurious low-depth noise (< min_germline_depth_ratio) is rejected
    call_noise = classify_pde_junction(p_noise_stem, graph, ledger)
    assert call_noise.is_pde_junction is False

    # Negative control 4: Heterozygous split where both branches have ~0.50x core depth (exceeding max_germline_ratio 0.33)
    graph.add_segment("het_stem", length=100_000, depth=60.0)
    graph.add_segment("het_arm_a", length=40_000, depth=30.0)
    graph.add_segment("het_arm_b", length=40_000, depth=30.0)

    p_het_stem = SegmentPort("het_stem", Strand.FORWARD)
    p_het_a = SegmentPort("het_arm_a", Strand.FORWARD)
    p_het_b = SegmentPort("het_arm_b", Strand.FORWARD)
    graph.add_edge(p_het_stem, p_het_a, cigar="100M")
    graph.add_edge(p_het_stem, p_het_b, cigar="100M")

    ev_het_stem = _make_evidence("het_stem", "3p", depth=60.0, degree=2, is_tip=False, usco_hit_count=8)
    ev_het_a = _make_evidence("het_arm_a", "3p", depth=30.0, degree=0, is_tip=True)
    ev_het_b = _make_evidence("het_arm_b", "3p", depth=30.0, degree=0, is_tip=True)

    ledger[p_het_stem.to_symbol()] = ev_het_stem
    ledger[p_het_a.to_symbol()] = ev_het_a
    ledger[p_het_b.to_symbol()] = ev_het_b

    call_het = classify_pde_junction(p_het_stem, graph, ledger)
    assert call_het.is_pde_junction is False, "Heterozygous split without stoichiometric GRS branch must be rejected"


# ==============================================================================
# BIO-004: Repeat Representation Preservation for Variable Satellites
# ==============================================================================


@pytest.mark.bio
def test_bio_004_variable_tandem_repeats_preserve_representative_sequences(tmp_path: Path):
    """BIO-004: Variable tandem repeats preserve representative sequences.

    Where variable-length tandem satellite arrays exhibit read disagreements across alleles,
    representative repeat sequences must be retained rather than collapsed into generic 100-bp N gaps.
    """
    graph = BidirectedGraph()
    monomer = "ACGTTAGGCT"
    representative_repeat = monomer * 25  # 250 bp of authentic tandem repeat

    graph.add_segment("stem", length=80_000, depth=64.0, sequence="A" * 80_000)
    graph.add_segment("som_cap", length=20_000, depth=45.0, sequence="C" * 20_000)
    graph.add_segment("sat_repeat_entry", length=len(representative_repeat), depth=35.0, sequence=representative_repeat)
    graph.add_segment("sat_alt_copy", length=len(representative_repeat), depth=12.0, sequence=representative_repeat)
    graph.add_segment("sat_alt_copy2", length=len(representative_repeat), depth=10.0, sequence=representative_repeat)

    p_stem = SegmentPort("stem", Strand.FORWARD)
    p_som = SegmentPort("som_cap", Strand.FORWARD)
    p_sat = SegmentPort("sat_repeat_entry", Strand.FORWARD)
    p_alt = SegmentPort("sat_alt_copy", Strand.FORWARD)
    p_alt2 = SegmentPort("sat_alt_copy2", Strand.FORWARD)

    graph.add_edge(p_stem, p_som, cigar="100M")
    graph.add_edge(p_stem, p_sat, cigar="100M")
    # Multivalent branching into alternate repeat copies (degree >= 2 on repeat)
    graph.add_edge(p_sat, p_alt, cigar="50M")
    graph.add_edge(p_sat, p_alt2, cigar="50M")

    ev_stem = _make_evidence("stem", "3p", depth=64.0, degree=2, is_tip=False, usco_hit_count=8)
    ev_som = _make_evidence("som_cap", "3p", depth=45.0, degree=0, is_tip=True, has_terminal_telomere=True)
    ev_sat = _make_evidence(
        "sat_repeat_entry",
        "3p",
        depth=35.0,
        length=len(representative_repeat),
        degree=2,
        is_tip=False,
    )
    ev_alt = _make_evidence(
        "sat_alt_copy",
        "3p",
        depth=12.0,
        length=len(representative_repeat),
        degree=0,
        is_tip=True,
    )
    ev_alt2 = _make_evidence(
        "sat_alt_copy2",
        "3p",
        depth=10.0,
        length=len(representative_repeat),
        degree=0,
        is_tip=True,
    )

    ledger = {
        p_stem.to_symbol(): ev_stem,
        p_som.to_symbol(): ev_som,
        p_sat.to_symbol(): ev_sat,
        p_alt.to_symbol(): ev_alt,
        p_alt2.to_symbol(): ev_alt2,
    }

    call: PdeJunctionCall = classify_pde_junction(p_stem, graph, ledger)
    assert call.is_pde_junction is True
    assert call.somatic_port == p_som
    # Confirm repeat entry is retained (not discarded due to multivalency)
    assert call.is_repeat_entry is True
    assert call.germline_port == p_sat

    # Verify that the repeat sequence is intact and not collapsed into Ns
    seq = graph.get_sequence("sat_repeat_entry")
    assert seq == representative_repeat
    assert "N" not in seq
    assert len(seq) == len(monomer) * 25

    # Lossless serialization round-trip: sequence remains identical
    out_gfa = tmp_path / "sat_test.gfa"
    graph.write_gfa(out_gfa)
    reloaded = load_gfa(out_gfa, load_sequences=True)
    assert reloaded.get_sequence("sat_repeat_entry") == representative_repeat

    # Negative control: Tandem repeat that loops directly back onto stem is rejected
    graph_loop = BidirectedGraph()
    graph_loop.add_segment("loop_stem", length=80_000, depth=64.0)
    graph_loop.add_segment("loop_som_cap", length=20_000, depth=45.0)
    graph_loop.add_segment("loop_repeat", length=len(representative_repeat), depth=35.0, sequence=representative_repeat)

    p_loop_stem = SegmentPort("loop_stem", Strand.FORWARD)
    p_loop_som = SegmentPort("loop_som_cap", Strand.FORWARD)
    p_loop_rep = SegmentPort("loop_repeat", Strand.FORWARD)

    graph_loop.add_edge(p_loop_stem, p_loop_som, cigar="100M")
    graph_loop.add_edge(p_loop_stem, p_loop_rep, cigar="100M")
    # Repeat loops back to stem
    graph_loop.add_edge(p_loop_rep, p_loop_stem, cigar="50M")

    ev_loop_stem = _make_evidence("loop_stem", "3p", depth=64.0, degree=2, is_tip=False, usco_hit_count=8)
    ev_loop_som = _make_evidence("loop_som_cap", "3p", depth=45.0, degree=0, is_tip=True, has_terminal_telomere=True)
    ev_loop_rep = _make_evidence("loop_repeat", "3p", depth=35.0, degree=2, is_tip=False)

    ledger_loop = {
        p_loop_stem.to_symbol(): ev_loop_stem,
        p_loop_som.to_symbol(): ev_loop_som,
        p_loop_rep.to_symbol(): ev_loop_rep,
    }

    call_loop = classify_pde_junction(p_loop_stem, graph_loop, ledger_loop)
    assert call_loop.is_pde_junction is False, (
        "Repeat looping back to stem must not be accepted as an authentic GRS branch"
    )


# ==============================================================================
# BIO-005: Lossless Overlap and Provenance Invariance
# ==============================================================================


@pytest.mark.bio
def test_bio_005_lossless_graph_representation_preserves_dovetails_and_offsets(tmp_path: Path):
    """BIO-005: Lossless graph representation preserves dovetails, indels, and offsets.

    Exact CIGAR operations with indels, segment lengths (LN), read depths (rd/dp),
    overlap tags (L1/L2), and embedded read placements (A-lines) must be preserved without loss.
    """
    graph = BidirectedGraph()
    graph.add_segment("seg_alpha", length=50_000, depth=64.0, sequence="A" * 50_000, raw_tags=["LN:i:50000", "rd:i:64"])
    graph.add_segment("seg_beta", length=30_000, depth=32.0, sequence="T" * 30_000, raw_tags=["LN:i:30000", "rd:i:32"])

    p_alpha = SegmentPort("seg_alpha", Strand.FORWARD)
    p_beta = SegmentPort("seg_beta", Strand.FORWARD)

    # Complex dovetail with insertions and deletions: 45M2D50M1I10M
    cigar_complex = "45M2D50M1I10M"
    graph.add_edge(
        p_alpha,
        p_beta,
        cigar=cigar_complex,
        l1=50_000,
        l2=30_000,
        raw_tags=["L1:i:50000", "L2:i:30000"],
    )

    # Add Hifiasm A-line read placement with haplogroup tag
    a_line = "A\tseg_alpha\t1500\t+\tread_00123\t0\t15000\tHG:A:p"
    graph.a_lines.append(a_line)

    gfa_file = tmp_path / "lossless_test.gfa"
    graph.write_gfa(gfa_file)

    # Load and verify complete fidelity
    loaded_a_lines: list[str] = []
    reloaded_graph = load_gfa(gfa_file, load_sequences=True, a_lines_store=loaded_a_lines)

    assert "seg_alpha" in reloaded_graph
    assert "seg_beta" in reloaded_graph
    assert reloaded_graph.segments["seg_alpha"].length == 50_000
    assert reloaded_graph.segments["seg_alpha"].depth == 64.0

    outgoing = reloaded_graph.get_outgoing(p_alpha)
    assert len(outgoing) == 1
    edge = outgoing[0]
    assert edge.cigar == cigar_complex
    # CIGAR overlap computation must match
    assert cigar_overlap_length(edge.cigar) == cigar_overlap_length(cigar_complex)
    assert edge.l1 == 50_000
    assert edge.l2 == 30_000

    # Watson-Crick dual edge must be preserved and match
    dual_outgoing = reloaded_graph.get_outgoing(p_beta.invert())
    assert len(dual_outgoing) == 1
    assert dual_outgoing[0].cigar == cigar_complex

    # Read placement A-lines must be retained losslessly
    assert len(loaded_a_lines) == 1
    assert loaded_a_lines[0] == a_line
    assert len(reloaded_graph.a_lines) == 1
    assert reloaded_graph.a_lines[0] == a_line

    # Non-standard segment and edge tags (LN, rd, L1, L2) must be preserved losslessly
    assert reloaded_graph.raw_segment_tags.get("seg_alpha") == ["LN:i:50000", "rd:i:64"]
    assert reloaded_graph.raw_segment_tags.get("seg_beta") == ["LN:i:30000", "rd:i:32"]
    assert reloaded_graph.raw_edge_tags.get((p_alpha, p_beta)) == ["L1:i:50000", "L2:i:30000"]


# ==============================================================================
# BIO-006: Held-Out Generalization Guard
# ==============================================================================


@pytest.mark.bio
def test_bio_006_unitig_id_permutation_and_sorting_invariance():
    """BIO-006: Re-labeling/permuting unitig IDs does not change PDE junction calling or cap identification.

    Graph decision logic must depend purely on topology and stoichiometry, never on lexical
    sort order, unitig naming conventions (utg...), or species IDs.
    """
    # Baseline scenario
    g_base = BidirectedGraph()
    g_base.add_segment("stem_01", length=100_000, depth=64.0)
    g_base.add_segment("som_02", length=25_000, depth=48.0)
    g_base.add_segment("grs_03", length=50_000, depth=7.0)

    p_stem = SegmentPort("stem_01", Strand.FORWARD)
    p_som = SegmentPort("som_02", Strand.FORWARD)
    p_grs = SegmentPort("grs_03", Strand.FORWARD)

    g_base.add_edge(p_stem, p_som, cigar="100M")
    g_base.add_edge(p_stem, p_grs, cigar="100M")

    ev_stem = _make_evidence("stem_01", "3p", depth=64.0, degree=2, is_tip=False, usco_hit_count=8)
    ev_som = _make_evidence("som_02", "3p", depth=48.0, degree=0, is_tip=True, has_terminal_telomere=True)
    ev_grs = _make_evidence("grs_03", "3p", depth=7.0, degree=0, is_tip=True)

    ledger_base = {
        p_stem.to_symbol(): ev_stem,
        p_som.to_symbol(): ev_som,
        p_grs.to_symbol(): ev_grs,
    }

    call_base: PdeJunctionCall = classify_pde_junction(p_stem, g_base, ledger_base)
    assert call_base.is_pde_junction is True
    assert call_base.somatic_port == p_som
    assert call_base.germline_port == p_grs

    # Permutation 1: Invert lexical sort order completely
    mapping_perm1 = {
        "stem_01": "utg999999",  # Originally first, now lexicographically last
        "som_02": "utg000001",  # Originally second, now lexicographically first
        "grs_03": "utg555555",  # Intermediate
    }

    g_perm1 = BidirectedGraph()
    g_perm1.add_segment(mapping_perm1["grs_03"], length=50_000, depth=7.0)
    g_perm1.add_segment(mapping_perm1["som_02"], length=25_000, depth=48.0)
    g_perm1.add_segment(mapping_perm1["stem_01"], length=100_000, depth=64.0)

    p_stem_p1 = SegmentPort(mapping_perm1["stem_01"], Strand.FORWARD)
    p_som_p1 = SegmentPort(mapping_perm1["som_02"], Strand.FORWARD)
    p_grs_p1 = SegmentPort(mapping_perm1["grs_03"], Strand.FORWARD)

    # Add edges in reverse order
    g_perm1.add_edge(p_stem_p1, p_grs_p1, cigar="100M")
    g_perm1.add_edge(p_stem_p1, p_som_p1, cigar="100M")

    ev_stem_p1 = _make_evidence(mapping_perm1["stem_01"], "3p", depth=64.0, degree=2, is_tip=False, usco_hit_count=8)
    ev_som_p1 = _make_evidence(
        mapping_perm1["som_02"], "3p", depth=48.0, degree=0, is_tip=True, has_terminal_telomere=True
    )
    ev_grs_p1 = _make_evidence(mapping_perm1["grs_03"], "3p", depth=7.0, degree=0, is_tip=True)

    ledger_perm1 = {
        p_stem_p1.to_symbol(): ev_stem_p1,
        p_som_p1.to_symbol(): ev_som_p1,
        p_grs_p1.to_symbol(): ev_grs_p1,
    }

    call_perm1: PdeJunctionCall = classify_pde_junction(p_stem_p1, g_perm1, ledger_perm1)
    assert call_perm1.is_pde_junction is True
    assert call_perm1.somatic_port == p_som_p1
    assert call_perm1.germline_port == p_grs_p1
    assert abs(call_perm1.depth_ratio - call_base.depth_ratio) < 1e-6

    # Permutation 2: Arbitrary UUID-like names with mixed cases
    mapping_perm2 = {
        "stem_01": "chrX_core_chunk_Alpha",
        "som_02": "telomere_cap_Omega",
        "grs_03": "grs_eliminated_Mu",
    }

    g_perm2 = BidirectedGraph()
    for old_name, new_name in mapping_perm2.items():
        g_perm2.add_segment(new_name, length=g_base.segments[old_name].length, depth=g_base.segments[old_name].depth)

    p_stem_p2 = SegmentPort(mapping_perm2["stem_01"], Strand.FORWARD)
    p_som_p2 = SegmentPort(mapping_perm2["som_02"], Strand.FORWARD)
    p_grs_p2 = SegmentPort(mapping_perm2["grs_03"], Strand.FORWARD)

    g_perm2.add_edge(p_stem_p2, p_som_p2, cigar="100M")
    g_perm2.add_edge(p_stem_p2, p_grs_p2, cigar="100M")

    ev_stem_p2 = _make_evidence(mapping_perm2["stem_01"], "3p", depth=64.0, degree=2, is_tip=False, usco_hit_count=8)
    ev_som_p2 = _make_evidence(
        mapping_perm2["som_02"], "3p", depth=48.0, degree=0, is_tip=True, has_terminal_telomere=True
    )
    ev_grs_p2 = _make_evidence(mapping_perm2["grs_03"], "3p", depth=7.0, degree=0, is_tip=True)

    ledger_perm2 = {
        p_stem_p2.to_symbol(): ev_stem_p2,
        p_som_p2.to_symbol(): ev_som_p2,
        p_grs_p2.to_symbol(): ev_grs_p2,
    }

    call_perm2: PdeJunctionCall = classify_pde_junction(p_stem_p2, g_perm2, ledger_perm2)
    assert call_perm2.is_pde_junction is True
    assert call_perm2.somatic_port == p_som_p2
    assert call_perm2.germline_port == p_grs_p2
    assert abs(call_perm2.depth_ratio - call_base.depth_ratio) < 1e-6

    # Verify cap identification invariance: somatic cap ports are recognized as TELOMERIC_CAP regardless of naming/ordering
    assert classify_port(p_som, g_base, ev_som).terminal_feature == TerminalFeature.TELOMERIC_CAP
    assert classify_port(p_som_p1, g_perm1, ev_som_p1).terminal_feature == TerminalFeature.TELOMERIC_CAP
    assert classify_port(p_som_p2, g_perm2, ev_som_p2).terminal_feature == TerminalFeature.TELOMERIC_CAP

    # Verify non-cap classification invariance: stem and GRS ports remain NONE across permutations
    for p, g, ev in [(p_stem, g_base, ev_stem), (p_stem_p1, g_perm1, ev_stem_p1), (p_stem_p2, g_perm2, ev_stem_p2)]:
        assert classify_port(p, g, ev).terminal_feature == TerminalFeature.NONE
    for p, g, ev in [(p_grs, g_base, ev_grs), (p_grs_p1, g_perm1, ev_grs_p1), (p_grs_p2, g_perm2, ev_grs_p2)]:
        assert classify_port(p, g, ev).terminal_feature == TerminalFeature.NONE
