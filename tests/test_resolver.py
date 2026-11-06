"""Tests for ONT spanning-read resolution engine and satellite gap finishing."""

import gzip
from io import StringIO
from pathlib import Path

import pytest

from gfapdefinisher.models import SegmentPort, Strand
from gfapdefinisher.policy import DEFAULT_POLICY
from gfapdefinisher.quirks import (
    EndState,
    IntervalQuirk,
    emit_end_states_tsv,
    emit_quirks_tsv,
    parse_end_states_tsv,
    parse_quirks_tsv,
)
from gfapdefinisher.resolver import (
    OntBridgingEngine,
    PafHit,
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
