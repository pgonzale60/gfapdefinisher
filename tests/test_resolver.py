"""Tests for ONT spanning-read resolution engine and satellite gap finishing."""

import gzip
from io import StringIO
from pathlib import Path

import pytest

from gfapdefinisher.graph import BidirectedGraph, load_gfa
from gfapdefinisher.ledger import EndEvidence
from gfapdefinisher.models import SegmentPort, Strand
from gfapdefinisher.policy import DEFAULT_POLICY
from gfapdefinisher.resolutions import (
    EndState,
    EndStateRecord,
    IntervalQuirk,
    QuirkRecord,
    emit_end_states_tsv,
    emit_quirks_tsv,
    parse_end_states_tsv,
    parse_quirks_tsv,
)
from gfapdefinisher.resolver import (
    OntBridgingEngine,
    PafHit,
    load_reads_sequences,
    parse_paf_hits,
    reverse_complement,
)


def _make_paf_line(
    qname: str,
    qlen: int,
    qstart: int,
    qend: int,
    strand: str,
    tname: str,
    tlen: int,
    tstart: int,
    tend: int,
    alen: int = 1500,
    mapq: int = 60,
    nmatch: int = 1400,
) -> str:
    return (
        f"{qname}\t{qlen}\t{qstart}\t{qend}\t{strand}\t"
        f"{tname}\t{tlen}\t{tstart}\t{tend}\t{nmatch}\t{alen}\t{mapq}\ttp:A:P"
    )


def test_paf_hit_parsing_and_slack_properties():
    """Verify PAF format parsing and boundary slack computation."""
    line = _make_paf_line(
        qname="read001",
        qlen=20000,
        qstart=1000,
        qend=2500,
        strand="+",
        tname="IV_LEFT",
        tlen=5000,
        tstart=3400,
        tend=4900,
        alen=1500,
        mapq=60,
    )
    hit = PafHit.from_line(line)

    assert hit.qname == "read001"
    assert hit.qlen == 20000
    assert hit.qstart == 1000
    assert hit.qend == 2500
    assert hit.strand == "+"
    assert hit.tname == "IV_LEFT"
    assert hit.tlen == 5000
    assert hit.tstart == 3400
    assert hit.tend == 4900
    assert hit.alen == 1500
    assert hit.mapq == 60
    assert hit.tags.get("tp") == "tp:A:P"
    assert hit.left_slack == 100  # 5000 - 4900
    assert hit.right_slack == 3400

    # Malformed line rejection
    with pytest.raises(ValueError, match="expected >= 12 columns"):
        PafHit.from_line("short\tline")


def test_parse_paf_hits_stream(tmp_path: Path):
    """Verify streaming PAF hits across text and gzipped inputs."""
    lines = [
        "# Header comment line",
        "",
        _make_paf_line("r1", 10000, 100, 1600, "+", "L", 5000, 3500, 5000),
        _make_paf_line("r2", 12000, 200, 1700, "-", "R", 5000, 0, 1500),
    ]
    paf_content = "\n".join(lines)

    # Stream from StringIO
    hits = list(parse_paf_hits(StringIO(paf_content)))
    assert len(hits) == 2
    assert hits[0].qname == "r1"
    assert hits[1].qname == "r2"

    # Stream from gzipped file
    gz_file = tmp_path / "test.paf.gz"
    with gzip.open(gz_file, "wt", encoding="utf-8") as handle:
        handle.write(paf_content)

    gz_hits = list(parse_paf_hits(gz_file))
    assert len(gz_hits) == 2
    assert gz_hits[0].qname == "r1"
    assert gz_hits[1].qname == "r2"


def test_boundary_slack_filtering():
    """Verify boundary slack filtering on both Left and Right candidate ports."""
    engine = OntBridgingEngine()

    # Case 1: Left slack within tolerance (tend=4900, tlen=5000 -> slack=100 <= 1500)
    # Right slack within tolerance (tstart=50 -> slack=50 <= 1500)
    lines_valid = [
        _make_paf_line("read_ok", 15000, 1000, 2500, "+", "LEFT", 5000, 3400, 4900),
        _make_paf_line("read_ok", 15000, 4000, 5500, "+", "RIGHT", 5000, 50, 1550),
    ]
    spanning = engine.collect_spanning_reads(lines_valid, "LEFT", "RIGHT")
    assert len(spanning) == 1
    assert spanning[0].qname == "read_ok"
    assert spanning[0].gap_length == (4000 - 50) - (2500 + 100)

    # Case 2: Left slack exceeds tolerance (tend=3000, tlen=5000 -> slack=2000 > 1500)
    lines_left_excess = [
        _make_paf_line("read_bad", 15000, 1000, 2500, "+", "LEFT", 5000, 1500, 3000),
        _make_paf_line("read_bad", 15000, 4000, 5500, "+", "RIGHT", 5000, 50, 1550),
    ]
    assert len(engine.collect_spanning_reads(lines_left_excess, "LEFT", "RIGHT")) == 0

    # Case 3: Right slack exceeds tolerance (tstart=2000 > 1500)
    lines_right_excess = [
        _make_paf_line("read_bad", 15000, 1000, 2500, "+", "LEFT", 5000, 3500, 5000),
        _make_paf_line("read_bad", 15000, 4000, 5500, "+", "RIGHT", 5000, 2000, 3500),
    ]
    assert len(engine.collect_spanning_reads(lines_right_excess, "LEFT", "RIGHT")) == 0


def test_tandem_repeat_disambiguation_forward_and_reverse():
    """Verify innermost boundary hit detection disambiguates multivalent tandem satellite repeats."""
    engine = OntBridgingEngine()

    # Forward strand (+):
    # Left flank hit: 1000..3000 (t=3000..5000, slack=0)
    # Right flank multivalent hits:
    # Outer hit 1: 8000..10000 (t=0..2000, slack=0)
    # Innermost hit 2: 5000..7000 (t=0..2000, slack=0) -> closest to Left!
    fwd_lines = [
        _make_paf_line("read_fwd", 20000, 1000, 3000, "+", "LEFT", 5000, 3000, 5000),
        _make_paf_line("read_fwd", 20000, 8000, 10000, "+", "RIGHT", 5000, 0, 2000),
        _make_paf_line("read_fwd", 20000, 5000, 7000, "+", "RIGHT", 5000, 0, 2000),
    ]
    spanning_fwd = engine.collect_spanning_reads(fwd_lines, "LEFT", "RIGHT")
    assert len(spanning_fwd) == 1
    # Innermost hit 2 (5000..7000) chosen, gap is 5000 - 3000 = 2000 bp
    assert spanning_fwd[0].gap_length == 2000
    assert spanning_fwd[0].q_from == 3000
    assert spanning_fwd[0].q_to == 5000

    # Reverse strand (-):
    # In negative strand, Right precedes Left.
    # Right multivalent hits:
    # Outer hit 1: 1000..3000 (t=0..2000, slack=0)
    # Innermost hit 2: 4000..6000 (t=0..2000, slack=0) -> closest to Left!
    # Left hit: 8000..10000 (t=3000..5000, slack=0)
    rev_lines = [
        _make_paf_line("read_rev", 20000, 1000, 3000, "-", "RIGHT", 5000, 0, 2000),
        _make_paf_line("read_rev", 20000, 4000, 6000, "-", "RIGHT", 5000, 0, 2000),
        _make_paf_line("read_rev", 20000, 8000, 10000, "-", "LEFT", 5000, 3000, 5000),
    ]
    spanning_rev = engine.collect_spanning_reads(rev_lines, "LEFT", "RIGHT")
    assert len(spanning_rev) == 1
    # Innermost hit 2 (4000..6000) chosen, gap is 8000 - 6000 = 2000 bp
    assert spanning_rev[0].gap_length == 2000
    assert spanning_rev[0].q_from == 6000
    assert spanning_rev[0].q_to == 8000


def test_dovetail_collinearity_and_orientation_checks():
    """Verify rejection of strand mismatch and non-collinear overlaps."""
    engine = OntBridgingEngine()

    # Strand mismatch: Left is +, Right is -
    mismatch_lines = [
        _make_paf_line("r_mismatch", 15000, 1000, 2500, "+", "LEFT", 5000, 3500, 5000),
        _make_paf_line("r_mismatch", 15000, 4000, 5500, "-", "RIGHT", 5000, 0, 1500),
    ]
    assert len(engine.collect_spanning_reads(mismatch_lines, "LEFT", "RIGHT")) == 0

    # Inverted collinearity on + strand (Right before Left)
    inverted_lines = [
        _make_paf_line("r_inv", 15000, 5000, 6500, "+", "LEFT", 5000, 3500, 5000),
        _make_paf_line("r_inv", 15000, 1000, 2500, "+", "RIGHT", 5000, 0, 1500),
    ]
    assert len(engine.collect_spanning_reads(inverted_lines, "LEFT", "RIGHT")) == 0

    # Insufficient anchor length (< min_anchor_bp 1000)
    short_anchor_lines = [
        _make_paf_line("r_short", 15000, 1000, 1500, "+", "LEFT", 5000, 4500, 5000, alen=500),
        _make_paf_line("r_short", 15000, 3000, 4500, "+", "RIGHT", 5000, 0, 1500, alen=1500),
    ]
    assert len(engine.collect_spanning_reads(short_anchor_lines, "LEFT", "RIGHT")) == 0


def test_gap_distribution_statistics():
    """Verify empirical min, median, max, and IQR computation."""
    engine = OntBridgingEngine()

    # Empty inputs
    dist_empty = engine.compute_distribution([])
    assert dist_empty.count == 0
    assert dist_empty.min_bp == 0
    assert dist_empty.median_bp == 0
    assert dist_empty.max_bp == 0
    assert dist_empty.iqr_bp == 0

    # Single length
    dist_single = engine.compute_distribution([1500])
    assert dist_single.count == 1
    assert dist_single.min_bp == 1500
    assert dist_single.median_bp == 1500
    assert dist_single.max_bp == 1500
    assert dist_single.iqr_bp == 0

    # Multi-value distribution
    lengths = [1000, 1200, 1500, 1800, 2000]
    dist = engine.compute_distribution(lengths)
    assert dist.count == 5
    assert dist.min_bp == 1000
    assert dist.median_bp == 1500
    assert dist.max_bp == 2000
    # Q1=1100, Q3=1900 -> IQR = 800
    assert dist.iqr_bp > 0


def test_representative_sequence_selection_bio004():
    """Verify representative sequence chosen from read closest to empirical median gap (BIO-004)."""
    engine = OntBridgingEngine()

    # 3 candidate reads with distinct gap lengths: 1000 bp, 1500 bp, 2500 bp
    # Median gap length is 1500 bp
    paf_lines = [
        # read_a: gap = 1000 bp
        _make_paf_line("read_a", 15000, 1000, 2500, "+", "LEFT", 5000, 3500, 5000),
        _make_paf_line("read_a", 15000, 3500, 5000, "+", "RIGHT", 5000, 0, 1500),
        # read_b: gap = 1500 bp (matches median!)
        _make_paf_line("read_b", 15000, 1000, 2500, "+", "LEFT", 5000, 3500, 5000),
        _make_paf_line("read_b", 15000, 4000, 5500, "+", "RIGHT", 5000, 0, 1500),
        # read_c: gap = 2500 bp
        _make_paf_line("read_c", 15000, 1000, 2500, "+", "LEFT", 5000, 3500, 5000),
        _make_paf_line("read_c", 15000, 5000, 6500, "+", "RIGHT", 5000, 0, 1500),
    ]

    reads_dict = {
        "read_a": "A" * 15000,
        "read_b": "A" * 2500 + ("CAGT" * 375) + "T" * 9500,  # 1500 bp gap sequence
        "read_c": "G" * 15000,
    }

    p_left = SegmentPort("LEFT", Strand.FORWARD)
    p_right = SegmentPort("RIGHT", Strand.FORWARD)

    res = engine.resolve_satellite_gap(
        chrom="IV",
        side="3p",
        left_port=p_left,
        right_port=p_right,
        paf_source=paf_lines,
        reads_sequences=reads_dict,
        repeat_monomer_len_bp=4,
    )

    assert res.state == EndState.END_COMPLETE
    assert res.end_state_record.resolved is True
    assert res.quirk_record.representative_sequence_selected is True
    assert res.representative_read is not None
    assert res.representative_read.qname == "read_b"
    assert res.representative_sequence == "CAGT" * 375
    assert len(res.representative_sequence) == 1500
    assert res.quirk_record.repeat_monomer_len_bp == 4
    assert res.quirk_record.copy_number_range is not None


def test_representative_sequence_reverse_complement_selection():
    """Verify representative sequence extracted on reverse strand is reverse-complemented."""
    engine = OntBridgingEngine()

    gap_seq_fwd = "AACCGGTT" * 100  # 800 bp
    gap_seq_rc = reverse_complement(gap_seq_fwd)

    paf_lines = [
        # read_rev: gap = 800 bp on negative strand
        _make_paf_line("read_rev1", 15000, 1000, 2500, "-", "RIGHT", 5000, 0, 1500),
        _make_paf_line("read_rev1", 15000, 3300, 4800, "-", "LEFT", 5000, 3500, 5000),
        _make_paf_line("read_rev2", 15000, 1000, 2500, "-", "RIGHT", 5000, 0, 1500),
        _make_paf_line("read_rev2", 15000, 3300, 4800, "-", "LEFT", 5000, 3500, 5000),
    ]

    reads_dict = {
        "read_rev1": "N" * 2500 + gap_seq_rc + "N" * 11700,
        "read_rev2": "N" * 2500 + gap_seq_rc + "N" * 11700,
    }

    res = engine.resolve_satellite_gap(
        chrom="IV",
        side="3p",
        left_port="LEFT",
        right_port="RIGHT",
        paf_source=paf_lines,
        reads_sequences=reads_dict,
    )

    assert res.state == EndState.END_COMPLETE
    assert res.representative_sequence == gap_seq_fwd


def test_unbridged_orphan_cap_emission():
    """Verify unbridged gap resolves to EndState.END_GERMLINE_ORPHAN."""
    engine = OntBridgingEngine()

    # Flanks with zero spanning reads
    res = engine.resolve_satellite_gap(
        chrom="V",
        side="5p",
        left_port="V_LEFT",
        right_port="V_RIGHT",
        paf_source=[],
    )

    assert res.state == EndState.END_GERMLINE_ORPHAN
    assert res.end_state_record.resolved is False
    assert res.end_state_record.germline_port is None
    assert res.quirk_record.representative_sequence_selected is False
    assert res.quirk_record.agp_component_type == "U"
    assert res.quirk_record.spanning_read_count == 0


def test_tsv_ledger_emission_and_serde(tmp_path: Path):
    """Verify TSV SerDe round-trip of QuirkRecord and EndStateRecord instances."""
    engine = OntBridgingEngine()

    paf_lines = [
        _make_paf_line("r1", 10000, 1000, 2500, "+", "L", 5000, 3500, 5000),
        _make_paf_line("r1", 10000, 3500, 5000, "+", "R", 5000, 0, 1500),
        _make_paf_line("r2", 10000, 1000, 2500, "+", "L", 5000, 3500, 5000),
        _make_paf_line("r2", 10000, 3700, 5200, "+", "R", 5000, 0, 1500),
    ]

    res_bridged = engine.resolve_satellite_gap("IV", "3p", "L", "R", paf_source=paf_lines)
    res_orphan = engine.resolve_satellite_gap("V", "5p", "L", "R", paf_source=[])

    quirks = [res_bridged.quirk_record, res_orphan.quirk_record]
    end_states = [res_bridged.end_state_record, res_orphan.end_state_record]

    # Emit TSV files
    quirks_file = tmp_path / "quirks.tsv"
    end_states_file = tmp_path / "end_states.tsv"

    emit_quirks_tsv(quirks, quirks_file)
    emit_end_states_tsv(end_states, end_states_file)

    # Parse back and compare
    parsed_quirks = parse_quirks_tsv(quirks_file)
    parsed_end_states = parse_end_states_tsv(end_states_file)

    assert len(parsed_quirks) == 2
    assert parsed_quirks[0].feature_id == "sat_gap_variable:IV:3p:1"
    assert parsed_quirks[0].quirk_type == IntervalQuirk.SAT_GAP_VARIABLE
    assert parsed_quirks[0].spanning_read_count == 2
    assert parsed_quirks[0].representative_sequence_selected is True
    assert parsed_quirks[1].feature_id == "sat_gap_variable:V:5p:1"
    assert parsed_quirks[1].representative_sequence_selected is False

    assert len(parsed_end_states) == 2
    assert parsed_end_states[0].state == EndState.END_COMPLETE
    assert parsed_end_states[0].resolved is True
    assert parsed_end_states[1].state == EndState.END_GERMLINE_ORPHAN
    assert parsed_end_states[1].resolved is False


@pytest.mark.benchmark
def test_benchmark_df5120_resolution_iv_and_v():
    """Benchmark integration test resolving Chr IV 3p and Chr V 5p on benchmark dataset."""
    benchmark_paf = Path("data/benchmark/DF5120/ont_to_flanks.paf")
    if not benchmark_paf.exists():
        pytest.skip(f"Benchmark PAF file {benchmark_paf} not found.")

    engine = OntBridgingEngine(DEFAULT_POLICY)

    # 1. Resolve Chr IV 3p: IV_LEFT to IV_RIGHT
    res_iv = engine.resolve_satellite_gap(
        chrom="IV",
        side="3p",
        left_port="IV_LEFT",
        right_port="IV_RIGHT",
        paf_source=benchmark_paf,
    )

    assert res_iv.state == EndState.END_COMPLETE
    assert res_iv.end_state_record.resolved is True
    assert res_iv.quirk_record.spanning_read_count is not None
    assert res_iv.quirk_record.spanning_read_count >= 2
    assert res_iv.quirk_record.representative_sequence_selected is True
    assert res_iv.distribution is not None
    assert res_iv.distribution.median_bp > 0
    assert res_iv.representative_read is not None

    # 2. Resolve Chr V 5p: V_LEFT to V_RIGHT
    res_v = engine.resolve_satellite_gap(
        chrom="V",
        side="5p",
        left_port="V_LEFT",
        right_port="V_RIGHT",
        paf_source=benchmark_paf,
    )

    # V_RIGHT is unbridged due to short anchor length / lack of long bridging reads
    assert res_v.state == EndState.END_GERMLINE_ORPHAN
    assert res_v.end_state_record.resolved is False
    assert res_v.quirk_record.spanning_read_count == 0
    assert res_v.quirk_record.representative_sequence_selected is False


def test_collect_spanning_reads_generator_input():
    """Verify generator of PAF lines is fully consumed without dropping initial records."""
    engine = OntBridgingEngine()
    lines = [
        _make_paf_line("read_ok", 15000, 1000, 2500, "+", "LEFT", 5000, 3400, 4900),
        _make_paf_line("read_ok", 15000, 4000, 5500, "+", "RIGHT", 5000, 50, 1550),
    ]

    # Generator input
    gen = (line for line in lines)
    spanning = engine.collect_spanning_reads(gen, "LEFT", "RIGHT")
    assert len(spanning) == 1
    assert spanning[0].qname == "read_ok"


def test_mixed_strand_port_bridging():
    """Verify bridging when ports have mixed orientations (e.g. utg000126l- distal cap)."""
    engine = OntBridgingEngine()

    p_left = SegmentPort("LEFT", Strand.FORWARD)  # Exit at tlen (5000)
    p_right = SegmentPort("RIGHT", Strand.REVERSE)  # Entry at tlen (5000)

    # Forward walk: matches LEFT on +, RIGHT on -
    lines_fwd = [
        _make_paf_line("r_fwd", 15000, 1000, 2500, "+", "LEFT", 5000, 3500, 5000),
        _make_paf_line("r_fwd", 15000, 4000, 5500, "-", "RIGHT", 5000, 3500, 5000),
    ]
    spanning_fwd = engine.collect_spanning_reads(lines_fwd, p_left, p_right)
    assert len(spanning_fwd) == 1
    assert spanning_fwd[0].strand == "+"
    assert spanning_fwd[0].gap_length == 1500

    # Reverse walk: matches RIGHT on +, LEFT on -
    lines_rev = [
        _make_paf_line("r_rev", 15000, 1000, 2500, "+", "RIGHT", 5000, 3500, 5000),
        _make_paf_line("r_rev", 15000, 4000, 5500, "-", "LEFT", 5000, 3500, 5000),
    ]
    spanning_rev = engine.collect_spanning_reads(lines_rev, p_left, p_right)
    assert len(spanning_rev) == 1
    assert spanning_rev[0].strand == "-"
    assert spanning_rev[0].gap_length == 1500


def test_resolve_satellite_gap_sequence_population_and_index():
    """Verify representative_read.sequence is populated and custom gap_index is recorded."""
    engine = OntBridgingEngine()

    paf_lines = [
        _make_paf_line("r1", 10000, 1000, 2500, "+", "L", 5000, 3500, 5000),
        _make_paf_line("r1", 10000, 3500, 5000, "+", "R", 5000, 0, 1500),
        _make_paf_line("r2", 10000, 1000, 2500, "+", "L", 5000, 3500, 5000),
        _make_paf_line("r2", 10000, 3500, 5000, "+", "R", 5000, 0, 1500),
    ]
    reads_seqs = {
        "r1": "A" * 2500 + "G" * 1000 + "T" * 6500,
        "r2": "A" * 2500 + "G" * 1000 + "T" * 6500,
    }

    # Pass pre-collected spanning reads without sequence, then pass sequences to resolve_satellite_gap
    spanning_without_seq = engine.collect_spanning_reads(paf_lines, "L", "R")
    assert spanning_without_seq[0].sequence is None

    res = engine.resolve_satellite_gap(
        chrom="IV",
        side="3p",
        left_port="L",
        right_port="R",
        gap_index=2,
        spanning_reads=spanning_without_seq,
        reads_sequences=reads_seqs,
        somatic_depth_ratio=0.98,
    )

    assert res.state == EndState.END_COMPLETE
    assert res.quirk_record.feature_id == "sat_gap_variable:IV:3p:2"
    assert res.representative_sequence == "G" * 1000
    assert res.representative_read is not None
    assert res.representative_read.sequence == "G" * 1000
    assert res.end_state_record.somatic_depth_ratio == 0.98


def test_resolve_somatic_dropout_spanning_rescue(tmp_path: Path):
    """Verify ONT-mediated rescue transitions END_SOMATIC_MISSING to END_COMPLETE."""
    # Synthetic assembly graph: Chr V 3p architecture
    # utg000004l+ stem bifurcates into utg000179l- (dropout stub) and utg000139l+ (germline)
    graph = BidirectedGraph()
    graph.add_segment("utg000004l", length=6949969, depth=64.0)
    graph.add_segment("utg000179l", length=19005, depth=1.0)
    graph.add_segment("utg000139l", length=53655, depth=8.0)

    p_stem = SegmentPort("utg000004l", Strand.FORWARD)
    p_som = SegmentPort("utg000179l", Strand.REVERSE)
    p_grs = SegmentPort("utg000139l", Strand.FORWARD)

    graph.add_edge(p_stem, p_som, cigar="100M")
    graph.add_edge(p_stem, p_grs, cigar="100M")

    # 3 candidate spanning reads with distinct gap lengths: 1000, 1500, 2000 bp
    # Median gap length is 1500 bp (read_b)
    paf_lines = [
        # read_a: gap = 1000 bp
        _make_paf_line("read_a", 25000, 1000, 3000, "+", "utg000004l", 6949969, 6947969, 6949969, alen=2000),
        _make_paf_line("read_a", 25000, 4000, 6000, "+", "utg000139l", 53655, 0, 2000, alen=2000),
        # read_b: gap = 1500 bp (median candidate)
        _make_paf_line("read_b", 25000, 1000, 3000, "+", "utg000004l", 6949969, 6947969, 6949969, alen=2000),
        _make_paf_line("read_b", 25000, 4500, 6500, "+", "utg000139l", 53655, 0, 2000, alen=2000),
        # read_c: gap = 2000 bp
        _make_paf_line("read_c", 25000, 1000, 3000, "+", "utg000004l", 6949969, 6947969, 6949969, alen=2000),
        _make_paf_line("read_c", 25000, 5000, 7000, "+", "utg000139l", 53655, 0, 2000, alen=2000),
    ]

    paf_file = tmp_path / "spanning_dropout.paf"
    paf_file.write_text("\n".join(paf_lines) + "\n", encoding="utf-8")

    rep_expected = "CATG" * 375  # 1500 bp
    reads_fasta = {
        "read_a": "A" * 3000 + ("T" * 1000) + "G" * 21000,
        "read_b": "A" * 3000 + rep_expected + "G" * 20500,
        "read_c": "A" * 3000 + ("C" * 2000) + "G" * 20000,
    }

    engine = OntBridgingEngine()
    quirk_rec, end_rec, rep_seq = engine.resolve_somatic_dropout(
        chrom="V",
        side="3p",
        stem_port=p_stem,
        germline_port=p_grs,
        somatic_port=p_som,
        paf_path=paf_file,
        reads_fasta=reads_fasta,
        core_coverage_bp_median=64.0,
        somatic_depth_ratio=0.015,
        germline_depth_ratio=0.125,
    )

    # 1. End state transitions from END_SOMATIC_MISSING to END_COMPLETE with resolved=True
    assert end_rec.state == EndState.END_COMPLETE
    assert end_rec.resolved is True
    assert end_rec.notes == "ONT-derived sequence resolved germline telomere"
    assert end_rec.somatic_port == p_som.to_symbol()
    assert end_rec.germline_port == p_grs.to_symbol()
    assert end_rec.shared_port is None
    assert end_rec.core_coverage_bp_median == 64.0

    # 2. Representative sequence selected matching median candidate read
    assert rep_seq == rep_expected
    assert len(rep_seq) == 1500

    # 3. Quirk record emitted documenting LOW_SUPPORT_STUB with representative_sequence_selected=True
    assert quirk_rec.quirk_type == IntervalQuirk.LOW_SUPPORT_STUB
    assert quirk_rec.feature_id == "low_support_stub:V:3p:1"
    assert quirk_rec.unitig_ids == ("utg000179l",)
    assert quirk_rec.representative_sequence_selected is True
    assert quirk_rec.spanning_read_count == 3
    assert quirk_rec.length_median_bp == 1500
    assert quirk_rec.left_anchor_port == p_stem.to_symbol()
    assert quirk_rec.right_anchor_port == p_grs.to_symbol()
    assert quirk_rec.agp_component_type == "W"
    assert quirk_rec.agp_linkage_evidence == "ont-spanning"
    assert "unpolished" in quirk_rec.notes
    assert "ONT-derived" in quirk_rec.notes


def test_resolve_somatic_dropout_unbridged_retention(tmp_path: Path):
    """Verify unbridged somatic dropout retains END_SOMATIC_MISSING with resolved=False."""
    graph = BidirectedGraph()
    graph.add_segment("utg000004l", length=6949969, depth=64.0)
    graph.add_segment("utg000179l", length=19005, depth=1.0)
    graph.add_segment("utg000139l", length=53655, depth=8.0)

    p_stem = SegmentPort("utg000004l", Strand.FORWARD)
    p_som = SegmentPort("utg000179l", Strand.REVERSE)
    p_grs = SegmentPort("utg000139l", Strand.FORWARD)

    # Empty PAF file (no spanning reads)
    empty_paf = tmp_path / "empty.paf"
    empty_paf.write_text("", encoding="utf-8")

    engine = OntBridgingEngine()
    quirk_rec, end_rec, rep_seq = engine.resolve_somatic_dropout(
        chrom="V",
        side="3p",
        stem_port=p_stem,
        germline_port=p_grs,
        somatic_port=p_som,
        paf_path=empty_paf,
        core_coverage_bp_median=64.0,
        somatic_depth_ratio=0.015,
        germline_depth_ratio=0.125,
    )

    # 1. Retains END_SOMATIC_MISSING with resolved=False
    assert end_rec.state == EndState.END_SOMATIC_MISSING
    assert end_rec.resolved is False
    assert end_rec.somatic_port is None
    assert end_rec.germline_port == p_grs.to_symbol()
    assert "Unrescued somatic dropout" in end_rec.notes

    # 2. No sequence extracted
    assert rep_seq is None

    # 3. Quirk record documents unrescued stub with representative_sequence_selected=False
    assert quirk_rec.quirk_type == IntervalQuirk.LOW_SUPPORT_STUB
    assert quirk_rec.feature_id == "low_support_stub:V:3p:1"
    assert quirk_rec.unitig_ids == ("utg000179l",)
    assert quirk_rec.representative_sequence_selected is False
    assert quirk_rec.spanning_read_count == 0
    assert quirk_rec.agp_component_type == "U"
    assert quirk_rec.agp_linkage_evidence == "unbridged"


def test_resolve_somatic_dropout_reverse_strand_extraction(tmp_path: Path):
    """Verify representative sequence extracted on reverse strand is reverse-complemented."""
    p_stem = SegmentPort("utg_stem", Strand.FORWARD)
    p_som = SegmentPort("utg_som", Strand.REVERSE)
    p_grs = SegmentPort("utg_grs", Strand.FORWARD)

    expected_gap_fwd = "AAGCTTCC" * 125  # 1000 bp
    expected_gap_rc = reverse_complement(expected_gap_fwd)

    # Negative strand walk: Right hit precedes Left hit
    paf_lines = [
        _make_paf_line("rev_read", 20000, 2000, 4000, "-", "utg_grs", 50000, 0, 2000, alen=2000),
        _make_paf_line("rev_read", 20000, 5000, 7000, "-", "utg_stem", 50000, 48000, 50000, alen=2000),
    ]

    paf_file = tmp_path / "rev.paf"
    paf_file.write_text("\n".join(paf_lines) + "\n", encoding="utf-8")

    reads_fasta = {
        "rev_read": "N" * 4000 + expected_gap_rc + "N" * 15000,
    }

    engine = OntBridgingEngine()
    quirk_rec, end_rec, rep_seq = engine.resolve_somatic_dropout(
        chrom="IV",
        side="5p",
        stem_port=p_stem,
        germline_port=p_grs,
        somatic_port=p_som,
        paf_path=paf_file,
        reads_fasta=reads_fasta,
    )

    assert end_rec.state == EndState.END_COMPLETE
    assert end_rec.resolved is True
    assert quirk_rec.representative_sequence_selected is True
    assert rep_seq == expected_gap_fwd
    assert len(rep_seq) == 1000


def test_resolve_somatic_dropout_nonexistent_paf():
    """Verify nonexistent PAF path cleanly retains END_SOMATIC_MISSING without crashing."""
    p_stem = SegmentPort("utg_stem", Strand.FORWARD)
    p_som = SegmentPort("utg_som", Strand.REVERSE)
    p_grs = SegmentPort("utg_grs", Strand.FORWARD)

    engine = OntBridgingEngine()
    quirk_rec, end_rec, rep_seq = engine.resolve_somatic_dropout(
        chrom="X",
        side="3p",
        stem_port=p_stem,
        germline_port=p_grs,
        somatic_port=p_som,
        paf_path="nonexistent_file_path_12345.paf",
    )

    assert end_rec.state == EndState.END_SOMATIC_MISSING
    assert end_rec.resolved is False
    assert rep_seq is None
    assert quirk_rec.representative_sequence_selected is False
    assert quirk_rec.spanning_read_count == 0


def test_resolve_somatic_dropout_distal_germline_cap_rescue(tmp_path: Path):
    """Verify ONT reads spanning from stem to distal germline telomeric cap rescue somatic dropout."""
    graph = BidirectedGraph()
    graph.add_segment("utg_stem", length=100000, depth=60.0)
    graph.add_segment("utg_som", length=15000, depth=1.0)
    graph.add_segment("utg_germ_mid", length=20000, depth=8.0)
    graph.add_segment("utg_germ_cap", length=30000, depth=8.0)

    p_stem = SegmentPort("utg_stem", Strand.FORWARD)
    p_som = SegmentPort("utg_som", Strand.REVERSE)
    p_germ_mid = SegmentPort("utg_germ_mid", Strand.FORWARD)
    p_germ_cap = SegmentPort("utg_germ_cap", Strand.FORWARD)

    # Stem bifurcates into som and germ_mid
    graph.add_edge(p_stem, p_som, cigar="100M")
    graph.add_edge(p_stem, p_germ_mid, cigar="100M")
    # germ_mid extends into germ_cap
    graph.add_edge(p_germ_mid.invert(), p_germ_cap, cigar="100M")

    ev_cap = EndEvidence(
        segment_name="utg_germ_cap",
        end="3p",
        port=p_germ_cap,
        length=30000,
        depth=8.0,
        degree=0,
        is_tip=True,
        outgoing_neighbors=(),
        has_terminal_telomere=True,
        telomere_canonical_count=150,
    )
    ledger = {p_germ_cap.to_symbol(): ev_cap}

    # PAF contains spanning reads from utg_stem to distal cap utg_germ_cap (not utg_germ_mid)
    paf_lines = [
        _make_paf_line("distal_read1", 40000, 1000, 3000, "+", "utg_stem", 100000, 98000, 100000, alen=2000),
        _make_paf_line("distal_read1", 40000, 4200, 6200, "+", "utg_germ_cap", 30000, 0, 2000, alen=2000),
        _make_paf_line("distal_read2", 40000, 1000, 3000, "+", "utg_stem", 100000, 98000, 100000, alen=2000),
        _make_paf_line("distal_read2", 40000, 4200, 6200, "+", "utg_germ_cap", 30000, 0, 2000, alen=2000),
    ]
    paf_file = tmp_path / "distal_cap.paf"
    paf_file.write_text("\n".join(paf_lines) + "\n", encoding="utf-8")

    reads_fasta = {
        "distal_read1": "A" * 3000 + "G" * 1200 + "T" * 35800,
        "distal_read2": "A" * 3000 + "G" * 1200 + "T" * 35800,
    }

    engine = OntBridgingEngine(graph=graph, ledger=ledger)
    quirk_rec, end_rec, rep_seq = engine.resolve_somatic_dropout(
        chrom="I",
        side="3p",
        stem_port=p_stem,
        germline_port=p_germ_mid,
        somatic_port=p_som,
        paf_path=paf_file,
        reads_fasta=reads_fasta,
        core_coverage_bp_median=60.0,
        somatic_depth_ratio=0.016,
        germline_depth_ratio=0.133,
    )

    assert end_rec.state == EndState.END_COMPLETE
    assert end_rec.resolved is True
    assert end_rec.notes == "ONT-derived sequence resolved germline telomere"
    assert end_rec.germline_port == p_germ_cap.to_symbol()
    assert quirk_rec.representative_sequence_selected is True
    assert quirk_rec.agp_linkage_evidence == "ont-spanning"
    assert "unpolished" in quirk_rec.notes
    assert "ONT-derived" in quirk_rec.notes
    assert rep_seq == "G" * 1200


def test_load_reads_sequences_selective_streaming(tmp_path: Path):
    """Verify load_reads_sequences extracts only specified target_qnames without buffering all reads."""
    # Test FASTA
    fa_path = tmp_path / "test_reads.fa"
    fa_content = (
        "\n".join(
            [
                ">read_A first read",
                "AAAAA",
                ">read_B second read",
                "CCCCC",
                ">read_C third read",
                "GGGGG",
                ">read_D fourth read",
                "TTTTT",
            ]
        )
        + "\n"
    )
    fa_path.write_text(fa_content, encoding="utf-8")

    # Selective extraction
    sub_fa = load_reads_sequences(fa_path, target_qnames={"read_B", "read_D"})
    assert sub_fa == {"read_B": "CCCCC", "read_D": "TTTTT"}

    # Target not present
    empty_fa = load_reads_sequences(fa_path, target_qnames={"nonexistent"})
    assert empty_fa == {}

    # Empty target set
    zero_fa = load_reads_sequences(fa_path, target_qnames=set())
    assert zero_fa == {}

    # None target extracts all
    all_fa = load_reads_sequences(fa_path)
    assert len(all_fa) == 4

    # Test FASTQ
    fq_path = tmp_path / "test_reads.fq"
    fq_content = (
        "\n".join(
            [
                "@fq_1 read one",
                "ACGTACGT",
                "+",
                "IIIIIIII",
                "@fq_2 read two",
                "TGCATGCA",
                "+",
                "JJJJJJJJ",
                "@fq_3 read three",
                "GGGGGGGG",
                "+",
                "KKKKKKKK",
            ]
        )
        + "\n"
    )
    fq_path.write_text(fq_content, encoding="utf-8")

    sub_fq = load_reads_sequences(fq_path, target_qnames=["fq_2"])
    assert sub_fq == {"fq_2": "TGCATGCA"}


def test_collect_spanning_read_ids(tmp_path: Path):
    """Verify collect_spanning_read_ids returns unique read IDs bridging candidate ports."""
    paf_lines = [
        _make_paf_line("span_1", 20000, 1000, 2500, "+", "utg_left", 30000, 28500, 30000),
        _make_paf_line("span_1", 20000, 4000, 5500, "+", "utg_right", 10000, 0, 1500),
        _make_paf_line("span_2", 20000, 1000, 2500, "+", "utg_left", 30000, 28500, 30000),
        _make_paf_line("span_2", 20000, 4000, 5500, "+", "utg_right", 10000, 0, 1500),
        _make_paf_line("unrelated_1", 20000, 1000, 2500, "+", "utg_left", 30000, 0, 1500),
    ]
    engine = OntBridgingEngine()
    target_ids = engine.collect_spanning_read_ids(paf_lines)
    assert target_ids == {"span_1", "span_2"}


def test_benchmark_df5120_identify_orphan_caps():
    """Verify that identify_orphan_caps on DF5120 benchmark dataset returns the 10 isolated caps."""
    gfa_path = Path("data/benchmark/DF5120/asm.p_utg.gfa.gz")
    if not gfa_path.exists():
        gfa_path = Path("data/benchmark/DF5120/asm.p_utg.gfa")
    ledger_path = Path("data/benchmark/DF5120/df5120_end_evidence_ledger.tsv")

    if not gfa_path.exists() or not ledger_path.exists():
        pytest.skip("Benchmark DF5120 dataset not present in data/benchmark/DF5120.")

    graph = load_gfa(gfa_path)
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

    engine = OntBridgingEngine(graph=graph, ledger=ledger)
    orphan_caps = engine.identify_orphan_caps()

    expected_caps = {
        "utg000053l",
        "utg000056l",
        "utg000098l",
        "utg000126l",
        "utg000150l",
        "utg000153l",
        "utg000160l",
        "utg000197l",
        "utg000213l",
        "utg000217l",
    }
    assert set(orphan_caps.keys()) == expected_caps
    assert "utg000009l" not in orphan_caps


def test_identify_and_resolve_gaps_prevents_overwrite_on_shared_target():
    """Verify identify_and_resolve_gaps does not overwrite an already-matched EndStateRecord.

    When multiple spanning pairs share a target name (e.g. chr_grs), the first resolved
    end must be retained and not overwritten by subsequent pairs sharing the target.
    """
    graph = BidirectedGraph()
    graph.add_segment("chr_stem", 50000)
    graph.add_segment("chr_som", 25000)
    graph.add_segment("chr_grs", 30000)
    graph.add_segment("orphan_cap1", 10000)
    graph.add_segment("orphan_cap2", 10000)

    ev1_plus = EndEvidence(
        segment_name="orphan_cap1",
        end="3p",
        port=SegmentPort("orphan_cap1", Strand.PLUS),
        length=10000,
        depth=5.0,
        degree=0,
        is_tip=True,
        outgoing_neighbors=(),
        has_terminal_telomere=True,
        usco_hit_count=0,
    )
    ev1_minus = EndEvidence(
        segment_name="orphan_cap1",
        end="5p",
        port=SegmentPort("orphan_cap1", Strand.MINUS),
        length=10000,
        depth=5.0,
        degree=0,
        is_tip=True,
        outgoing_neighbors=(),
        has_terminal_telomere=False,
        usco_hit_count=0,
    )
    ev2_plus = EndEvidence(
        segment_name="orphan_cap2",
        end="3p",
        port=SegmentPort("orphan_cap2", Strand.PLUS),
        length=10000,
        depth=5.0,
        degree=0,
        is_tip=True,
        outgoing_neighbors=(),
        has_terminal_telomere=True,
        usco_hit_count=0,
    )
    ev2_minus = EndEvidence(
        segment_name="orphan_cap2",
        end="5p",
        port=SegmentPort("orphan_cap2", Strand.MINUS),
        length=10000,
        depth=5.0,
        degree=0,
        is_tip=True,
        outgoing_neighbors=(),
        has_terminal_telomere=False,
        usco_hit_count=0,
    )

    ledger = {
        "orphan_cap1+": ev1_plus,
        "orphan_cap1-": ev1_minus,
        "orphan_cap2+": ev2_plus,
        "orphan_cap2-": ev2_minus,
    }

    engine = OntBridgingEngine(graph=graph, ledger=ledger)
    orphan_caps = engine.identify_orphan_caps()
    assert "orphan_cap1" in orphan_caps
    assert "orphan_cap2" in orphan_caps

    initial_end = EndStateRecord(
        feature_id="end_chr_stem_3p",
        chrom="chr_stem",
        side="3p",
        state=EndState.END_SOMATIC_MISSING,
        somatic_port="chr_som+",
        germline_port="chr_grs+",
        resolved=False,
    )
    initial_quirk = QuirkRecord(
        feature_id="q_chr_stem_3p",
        chrom="chr_stem",
        side="3p",
        quirk_type=IntervalQuirk.SAT_GAP_VARIABLE,
    )

    # Pair 1: chr_grs to orphan_cap1 (3 spanning reads)
    # Pair 2: chr_grs to orphan_cap2 (3 spanning reads, sharing chr_grs target)
    paf_lines = [
        _make_paf_line("r1", 20000, 1000, 2500, "+", "chr_grs", 30000, 28500, 30000),
        _make_paf_line("r1", 20000, 4000, 5500, "+", "orphan_cap1", 10000, 0, 1500),
        _make_paf_line("r2", 20000, 1000, 2500, "+", "chr_grs", 30000, 28500, 30000),
        _make_paf_line("r2", 20000, 4000, 5500, "+", "orphan_cap1", 10000, 0, 1500),
        _make_paf_line("r3", 20000, 1000, 2500, "+", "chr_grs", 30000, 28500, 30000),
        _make_paf_line("r3", 20000, 4000, 5500, "+", "orphan_cap1", 10000, 0, 1500),
        _make_paf_line("r4", 20000, 1000, 2500, "+", "chr_grs", 30000, 28500, 30000),
        _make_paf_line("r4", 20000, 4000, 5500, "+", "orphan_cap2", 10000, 0, 1500),
        _make_paf_line("r5", 20000, 1000, 2500, "+", "chr_grs", 30000, 28500, 30000),
        _make_paf_line("r5", 20000, 4000, 5500, "+", "orphan_cap2", 10000, 0, 1500),
        _make_paf_line("r6", 20000, 1000, 2500, "+", "chr_grs", 30000, 28500, 30000),
        _make_paf_line("r6", 20000, 4000, 5500, "+", "orphan_cap2", 10000, 0, 1500),
    ]

    end_states, quirks = engine.identify_and_resolve_gaps(
        paf_source=paf_lines,
        end_states=[initial_end],
        quirks=[initial_quirk],
    )

    stem_ends = [e for e in end_states if e.chrom == "chr_stem" and e.side == "3p"]
    assert len(stem_ends) == 1
    resolved_stem = stem_ends[0]

    # Verify that chr_stem:3p resolved to orphan_cap1 and was NOT overwritten by orphan_cap2
    assert resolved_stem.state == EndState.END_COMPLETE
    assert resolved_stem.resolved is True
    assert resolved_stem.somatic_port == "chr_som+"
    assert resolved_stem.germline_port == "orphan_cap1+"


def test_identify_and_resolve_gaps_opposite_strand_cap():
    """Verify gap resolution when orphan cap terminal telomere is at 5' end (port cap-)."""
    graph = BidirectedGraph()
    graph.add_segment("chr_anchor", 50000)
    graph.add_segment("cap_rev", 15000)

    # cap_rev has terminal telomere at 5' end (port -), so port cap_rev- is telomeric
    ev_plus = EndEvidence(
        segment_name="cap_rev",
        end="3p",
        port=SegmentPort("cap_rev", Strand.PLUS),
        length=15000,
        depth=5.0,
        degree=0,
        is_tip=True,
        outgoing_neighbors=(),
        has_terminal_telomere=False,
        usco_hit_count=0,
    )
    ev_minus = EndEvidence(
        segment_name="cap_rev",
        end="5p",
        port=SegmentPort("cap_rev", Strand.MINUS),
        length=15000,
        depth=5.0,
        degree=0,
        is_tip=True,
        outgoing_neighbors=(),
        has_terminal_telomere=True,
        usco_hit_count=0,
    )

    ledger = {
        "cap_rev+": ev_plus,
        "cap_rev-": ev_minus,
    }

    engine = OntBridgingEngine(graph=graph, ledger=ledger)
    orphan_caps = engine.identify_orphan_caps()
    assert "cap_rev" in orphan_caps
    assert orphan_caps["cap_rev"] == ("cap_rev-", "cap_rev+")

    initial_end = EndStateRecord(
        feature_id="end_chr4_5p",
        chrom="IV",
        side="5p",
        state=EndState.END_SOMATIC_MISSING,
        somatic_port="chr4_som-",
        germline_port="chr_anchor+",
        resolved=False,
    )
    initial_quirk = QuirkRecord(
        feature_id="q_chr4_5p",
        chrom="IV",
        side="5p",
        quirk_type=IntervalQuirk.SAT_GAP_VARIABLE,
    )

    paf_lines = [
        _make_paf_line("r1", 20000, 1000, 2500, "+", "chr_anchor", 50000, 48500, 50000),
        _make_paf_line("r1", 20000, 4000, 5500, "-", "cap_rev", 15000, 13500, 15000),
        _make_paf_line("r2", 20000, 1000, 2500, "+", "chr_anchor", 50000, 48500, 50000),
        _make_paf_line("r2", 20000, 4000, 5500, "-", "cap_rev", 15000, 13500, 15000),
        _make_paf_line("r3", 20000, 1000, 2500, "+", "chr_anchor", 50000, 48500, 50000),
        _make_paf_line("r3", 20000, 4000, 5500, "-", "cap_rev", 15000, 13500, 15000),
    ]

    end_states, quirks = engine.identify_and_resolve_gaps(
        paf_source=paf_lines,
        end_states=[initial_end],
        quirks=[initial_quirk],
    )

    assert len(end_states) == 1
    res_end = end_states[0]
    assert res_end.state == EndState.END_COMPLETE
    assert res_end.resolved is True
    assert res_end.somatic_port == "chr4_som-"
    assert res_end.germline_port == "cap_rev-"

    assert len(quirks) == 1
    res_quirk = quirks[0]
    assert res_quirk.left_anchor_port == "chr_anchor+"
    assert res_quirk.right_anchor_port == "cap_rev-"
    assert res_quirk.spanning_read_count == 3


def test_identify_and_resolve_gaps_filters_unmatched_pairs_when_end_states_provided():
    """Verify that unassigned pairs in PAF do not create spurious ends when end_states is provided."""
    graph = BidirectedGraph()
    graph.add_segment("chr_anchor", 50000)
    graph.add_segment("unrelated_a", 20000)
    graph.add_segment("unrelated_b", 20000)

    engine = OntBridgingEngine(graph=graph, ledger={})

    initial_end = EndStateRecord(
        feature_id="end_chr1_3p",
        chrom="I",
        side="3p",
        state=EndState.END_SOMATIC_MISSING,
        somatic_port="chr1_som+",
        germline_port="chr_anchor+",
        resolved=False,
    )

    paf_lines = [
        _make_paf_line("r1", 20000, 1000, 2500, "+", "unrelated_a", 20000, 18500, 20000),
        _make_paf_line("r1", 20000, 4000, 5500, "+", "unrelated_b", 20000, 0, 1500),
        _make_paf_line("r2", 20000, 1000, 2500, "+", "unrelated_a", 20000, 18500, 20000),
        _make_paf_line("r2", 20000, 4000, 5500, "+", "unrelated_b", 20000, 0, 1500),
        _make_paf_line("r3", 20000, 1000, 2500, "+", "unrelated_a", 20000, 18500, 20000),
        _make_paf_line("r3", 20000, 4000, 5500, "+", "unrelated_b", 20000, 0, 1500),
    ]

    end_states, quirks = engine.identify_and_resolve_gaps(
        paf_source=paf_lines,
        end_states=[initial_end],
        quirks=[],
    )

    assert len(end_states) == 1
    assert end_states[0].chrom == "I"
    assert end_states[0].resolved is False
    assert len(quirks) == 0


def test_identify_and_resolve_gaps_two_hop_grs():
    """Verify 2-hop ONT GRS gap bridging: anchor -> mid_grs -> cap."""
    graph = BidirectedGraph()
    graph.add_segment("chr3_anchor", 50000)
    graph.add_segment("mid_grs", 25000)
    graph.add_segment("cap_3p", 10000)

    mid_plus = EndEvidence(
        segment_name="mid_grs",
        end="3p",
        port=SegmentPort("mid_grs", Strand.PLUS),
        length=25000,
        depth=5.0,
        degree=0,
        is_tip=True,
        outgoing_neighbors=(),
        has_terminal_telomere=False,
        usco_hit_count=0,
    )
    mid_minus = EndEvidence(
        segment_name="mid_grs",
        end="5p",
        port=SegmentPort("mid_grs", Strand.MINUS),
        length=25000,
        depth=5.0,
        degree=0,
        is_tip=True,
        outgoing_neighbors=(),
        has_terminal_telomere=False,
        usco_hit_count=0,
    )

    cap_plus = EndEvidence(
        segment_name="cap_3p",
        end="3p",
        port=SegmentPort("cap_3p", Strand.PLUS),
        length=10000,
        depth=5.0,
        degree=0,
        is_tip=True,
        outgoing_neighbors=(),
        has_terminal_telomere=True,
        usco_hit_count=0,
    )
    cap_minus = EndEvidence(
        segment_name="cap_3p",
        end="5p",
        port=SegmentPort("cap_3p", Strand.MINUS),
        length=10000,
        depth=5.0,
        degree=0,
        is_tip=True,
        outgoing_neighbors=(),
        has_terminal_telomere=False,
        usco_hit_count=0,
    )

    ledger = {
        "mid_grs+": mid_plus,
        "mid_grs-": mid_minus,
        "cap_3p+": cap_plus,
        "cap_3p-": cap_minus,
    }

    engine = OntBridgingEngine(graph=graph, ledger=ledger)
    orphan_caps = engine.identify_orphan_caps()
    assert "cap_3p" in orphan_caps

    unlinked_grs = engine.identify_unlinked_grs_unitigs()
    assert "mid_grs" in unlinked_grs

    initial_end = EndStateRecord(
        feature_id="end_chr3_3p",
        chrom="III",
        side="3p",
        state=EndState.END_SOMATIC_MISSING,
        somatic_port="chr3_som+",
        germline_port="chr3_anchor+",
        resolved=False,
    )
    initial_quirk = QuirkRecord(
        feature_id="q_chr3_3p",
        chrom="III",
        side="3p",
        quirk_type=IntervalQuirk.SAT_GAP_VARIABLE,
    )

    paf_lines = [
        # Hop 1 reads (r1, r2, r3)
        _make_paf_line("r1", 20000, 1000, 2500, "+", "chr3_anchor", 50000, 48500, 50000),
        _make_paf_line("r1", 20000, 4000, 5500, "+", "mid_grs", 25000, 0, 1500),
        _make_paf_line("r2", 20000, 1000, 2500, "+", "chr3_anchor", 50000, 48500, 50000),
        _make_paf_line("r2", 20000, 4000, 5500, "+", "mid_grs", 25000, 0, 1500),
        _make_paf_line("r3", 20000, 1000, 2500, "+", "chr3_anchor", 50000, 48500, 50000),
        _make_paf_line("r3", 20000, 4000, 5500, "+", "mid_grs", 25000, 0, 1500),
        # Hop 2 reads (r4, r5, r6)
        _make_paf_line("r4", 20000, 1000, 2500, "-", "mid_grs", 25000, 0, 1500),
        _make_paf_line("r4", 20000, 4000, 5500, "+", "cap_3p", 10000, 0, 1500),
        _make_paf_line("r5", 20000, 1000, 2500, "-", "mid_grs", 25000, 0, 1500),
        _make_paf_line("r5", 20000, 4000, 5500, "+", "cap_3p", 10000, 0, 1500),
        _make_paf_line("r6", 20000, 1000, 2500, "-", "mid_grs", 25000, 0, 1500),
        _make_paf_line("r6", 20000, 4000, 5500, "+", "cap_3p", 10000, 0, 1500),
    ]

    end_states, quirks = engine.identify_and_resolve_gaps(
        paf_source=paf_lines,
        end_states=[initial_end],
        quirks=[initial_quirk],
    )

    assert len(end_states) == 1
    res_end = end_states[0]
    assert res_end.state == EndState.END_COMPLETE
    assert res_end.resolved is True
    assert res_end.somatic_port == "chr3_som+"
    assert res_end.germline_port == "cap_3p+"

    assert len(quirks) == 2
    q1 = next(q for q in quirks if ":1" in q.feature_id)
    q2 = next(q for q in quirks if ":2" in q.feature_id)
    assert q1.left_anchor_port == "chr3_anchor+"
    assert q1.right_anchor_port == "mid_grs+"
    assert q2.left_anchor_port == "mid_grs-"
    assert q2.right_anchor_port == "cap_3p+"


def test_resolve_subtelomeric_satellite_caps_emits_hifi_spanning():
    """Verify resolve_subtelomeric_satellite_caps sets agp_linkage_evidence to 'hifi-spanning'."""
    graph = BidirectedGraph()
    rep_unit = "ACGTACGTACGTACG"  # 15-mer
    anchor_seq = ("C" * 1000) + (rep_unit * 200)
    cap_seq = (rep_unit * 10) + ("TTAGGC" * 35)

    graph.add_segment("anchor_sat", length=len(anchor_seq), sequence=anchor_seq, depth=40.0)
    graph.add_segment("cap_sat", length=len(cap_seq), sequence=cap_seq, depth=35.0)

    ledger = {
        "anchor_sat+": EndEvidence(
            segment_name="anchor_sat",
            end="3p",
            port=SegmentPort("anchor_sat", Strand.PLUS),
            length=len(anchor_seq),
            depth=40.0,
            degree=0,
            is_tip=True,
            outgoing_neighbors=(),
            has_terminal_telomere=False,
            usco_hit_count=0,
        ),
        "cap_sat+": EndEvidence(
            segment_name="cap_sat",
            end="3p",
            port=SegmentPort("cap_sat", Strand.PLUS),
            length=len(cap_seq),
            depth=35.0,
            degree=0,
            is_tip=True,
            outgoing_neighbors=(),
            has_terminal_telomere=True,
            usco_hit_count=0,
        ),
    }

    engine = OntBridgingEngine(graph=graph, ledger=ledger)

    er = EndStateRecord(
        feature_id="end_I_3p",
        chrom="I",
        side="3p",
        state=EndState.END_SOMATIC_MISSING,
        somatic_port="anchor_sat+",
        germline_port="anchor_sat+",
        resolved=False,
    )
    orphan_caps = {"cap_sat": ("cap_sat+", "3p")}

    end_states, quirks = engine.resolve_subtelomeric_satellite_caps(
        end_states=[er],
        quirks=[],
        orphan_caps=orphan_caps,
    )

    assert len(quirks) == 1
    q = quirks[0]
    assert q.agp_linkage_evidence == "hifi-spanning"
    assert q.agp_component_type == "U"
    assert q.length_median_bp == 100
    assert q.unitig_ids == ("anchor_sat", "cap_sat")
