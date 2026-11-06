"""Unit tests for native terminal telomere tract scanner."""

from gfapdefinisher.telomere import TerminalTelomereCall, TerminalTelomereResult, scan_terminal_tract


def test_exact_terminal_cap_at_offset_0_prefix():
    """Verify bona fide terminal cap starting at offset 0 on 5' prefix."""
    sequence = ("GCCTAA" * 40) + ("A" * 5000)
    call = scan_terminal_tract(sequence, motif="GCCTAA", is_suffix=False)

    assert isinstance(call, TerminalTelomereCall)
    assert isinstance(call, TerminalTelomereResult)
    assert call.offset_bp == 0
    assert call.repeat_count == 40
    assert call.tract_length_bp == 240
    assert call.is_terminal_cap is True
    assert call.motif == "GCCTAA"


def test_exact_terminal_cap_at_offset_0_suffix():
    """Verify bona fide terminal cap ending at physical 3' terminus on suffix."""
    sequence = ("A" * 5000) + ("TTAGGC" * 40)
    call = scan_terminal_tract(sequence, motif="TTAGGC", is_suffix=True)

    assert call.offset_bp == 0
    assert call.repeat_count == 40
    assert call.tract_length_bp == 240
    assert call.is_terminal_cap is True
    assert call.motif == "TTAGGC"


def test_interstitial_telomeric_sequence_classified_non_terminal_prefix():
    """Verify ITS tract with high repeat count located > max_terminal_offset_bp is classified as cap=False."""
    sequence = ("A" * 2000) + ("GCCTAA" * 40) + ("A" * 3000)
    call = scan_terminal_tract(sequence, motif="GCCTAA", is_suffix=False)

    assert call.offset_bp == 2000
    assert call.repeat_count == 40
    assert call.is_terminal_cap is False


def test_interstitial_telomeric_sequence_classified_non_terminal_suffix():
    """Verify ITS tract on 3' suffix located > max_terminal_offset_bp from 3' tip is cap=False."""
    sequence = ("A" * 3000) + ("TTAGGC" * 40) + ("A" * 2000)
    call = scan_terminal_tract(sequence, motif="TTAGGC", is_suffix=True)

    assert call.offset_bp == 2000
    assert call.repeat_count == 40
    assert call.is_terminal_cap is False


def test_terminal_offset_threshold_boundary():
    """Verify calibrated max_terminal_offset_bp = 25 gates cap status cleanly."""
    # Offset = 25 bp (tolerates <= 25 bp non-telomeric leader sequence)
    seq_25 = ("A" * 25) + ("GCCTAA" * 35) + ("A" * 1000)
    call_25 = scan_terminal_tract(seq_25, motif="GCCTAA", is_suffix=False)
    assert call_25.offset_bp == 25
    assert call_25.repeat_count == 35
    assert call_25.is_terminal_cap is True

    # Offset = 26 bp (exceeds 25 bp limit -> non-terminal / ITS)
    seq_26 = ("A" * 26) + ("GCCTAA" * 35) + ("A" * 1000)
    call_26 = scan_terminal_tract(seq_26, motif="GCCTAA", is_suffix=False)
    assert call_26.offset_bp == 26
    assert call_26.repeat_count == 35
    assert call_26.is_terminal_cap is False


def test_sub_threshold_repeat_count():
    """Verify repeat count below min_terminal_repeat_count (30) fails cap classification."""
    sequence = ("GCCTAA" * 29) + ("A" * 2000)
    call = scan_terminal_tract(sequence, motif="GCCTAA", is_suffix=False)

    assert call.offset_bp == 0
    assert call.repeat_count == 29
    assert call.is_terminal_cap is False


def test_exact_threshold_repeat_count():
    """Verify exact threshold 30 qualifies as a terminal cap."""
    sequence = ("GCCTAA" * 30) + ("A" * 2000)
    call = scan_terminal_tract(sequence, motif="GCCTAA", is_suffix=False)

    assert call.offset_bp == 0
    assert call.repeat_count == 30
    assert call.is_terminal_cap is True


def test_gapped_tract_chaining_within_max_internal_gap():
    """Verify non-repeat gaps <= max_internal_gap_bp (25) are chained into single tract."""
    # 20 repeats + 25 bp non-repeat gap + 20 repeats
    sequence = ("TTAGGC" * 20) + ("N" * 25) + ("TTAGGC" * 20) + ("A" * 2000)
    call = scan_terminal_tract(sequence, motif="TTAGGC", is_suffix=False)

    expected_len = (20 * 6) + 25 + (20 * 6)
    assert call.repeat_count == 40
    assert call.tract_length_bp == expected_len
    assert call.offset_bp == 0
    assert call.is_terminal_cap is True


def test_gapped_tract_broken_when_gap_exceeds_threshold():
    """Verify non-repeat gaps > max_internal_gap_bp (25) break into separate tracts."""
    # 20 repeats + 26 bp non-repeat gap + 20 repeats
    sequence = ("TTAGGC" * 20) + ("N" * 26) + ("TTAGGC" * 20) + ("A" * 2000)
    call = scan_terminal_tract(sequence, motif="TTAGGC", is_suffix=False)

    # Chaining breaks; both tracts have 20 repeats (< 30), so cap is False
    assert call.repeat_count == 20
    assert call.is_terminal_cap is False


def test_short_sequences_smaller_than_scan_window():
    """Verify scanning works seamlessly on sequences shorter than terminal_scan_bp."""
    seq_5p = "GCCTAA" * 35  # length 210 bp < 10000 bp
    call_5p = scan_terminal_tract(seq_5p, motif="GCCTAA", is_suffix=False)
    assert call_5p.offset_bp == 0
    assert call_5p.repeat_count == 35
    assert call_5p.is_terminal_cap is True

    seq_3p = "TTAGGC" * 35
    call_3p = scan_terminal_tract(seq_3p, motif="TTAGGC", is_suffix=True)
    assert call_3p.offset_bp == 0
    assert call_3p.repeat_count == 35
    assert call_3p.is_terminal_cap is True


def test_empty_and_noise_sequences():
    """Verify clean zero handling on empty or non-matching sequences."""
    call_empty = scan_terminal_tract("", motif="TTAGGC")
    assert call_empty.tract_length_bp == 0
    assert call_empty.repeat_count == 0
    assert call_empty.offset_bp == 0
    assert call_empty.is_terminal_cap is False

    call_noise = scan_terminal_tract("A" * 1000, motif="TTAGGC")
    assert call_noise.tract_length_bp == 0
    assert call_noise.repeat_count == 0
    assert call_noise.offset_bp == 0
    assert call_noise.is_terminal_cap is False


def test_case_insensitivity():
    """Verify sequence and motif matching is case-insensitive."""
    sequence = ("gcctaa" * 35) + ("a" * 500)
    call = scan_terminal_tract(sequence, motif="GCCTAA", is_suffix=False)
    assert call.repeat_count == 35
    assert call.is_terminal_cap is True


def test_custom_policy_parameter_passing():
    """Verify explicit custom scan parameters override defaults."""
    sequence = ("A" * 30) + ("GCCTAA" * 50) + ("A" * 500)
    # Default: max_terminal_offset_bp = 25 -> cap is False
    call_def = scan_terminal_tract(sequence, motif="GCCTAA", is_suffix=False)
    assert call_def.is_terminal_cap is False

    # Custom: max_terminal_offset_bp = 35 -> cap is True
    call_cust = scan_terminal_tract(
        sequence,
        motif="GCCTAA",
        terminal_scan_bp=5000,
        max_terminal_offset_bp=35,
        max_internal_gap_bp=20,
        min_repeat_count=40,
        is_suffix=False,
    )
    assert call_cust.is_terminal_cap is True
