"""Unit tests for CIGAR parsing, consumption calculation, and excision offsets."""

import pytest

from gfapdefinisher.cigar import (
    CigarOp,
    calculate_excision_offsets,
    cigar_overlap_length,
    cigar_query_consumption,
    cigar_target_consumption,
    format_cigar,
    parse_cigar,
)


def test_parse_cigar_simple():
    """Verify parsing a simple single-operation CIGAR string."""
    ops = parse_cigar("15643M")
    assert len(ops) == 1
    assert ops[0] == CigarOp(length=15643, op="M")
    assert format_cigar(ops) == "15643M"


def test_parse_cigar_indels():
    """Verify parsing CIGAR strings with insertion and deletion operations."""
    # Insertion
    ops_ins = parse_cigar("100M2I50M")
    assert ops_ins == [
        CigarOp(length=100, op="M"),
        CigarOp(length=2, op="I"),
        CigarOp(length=50, op="M"),
    ]
    assert format_cigar(ops_ins) == "100M2I50M"

    # Deletion
    ops_del = parse_cigar("200M3D10M")
    assert ops_del == [
        CigarOp(length=200, op="M"),
        CigarOp(length=3, op="D"),
        CigarOp(length=10, op="M"),
    ]
    assert format_cigar(ops_del) == "200M3D10M"


def test_parse_cigar_empty_and_special():
    """Verify edge cases: empty string, unmapped star, and complex CIGAR."""
    assert parse_cigar("") == []
    assert parse_cigar("*") == []

    # Sequence match (=) and mismatch (X)
    ops_diff = parse_cigar("50=2X48=")
    assert ops_diff == [
        CigarOp(length=50, op="="),
        CigarOp(length=2, op="X"),
        CigarOp(length=48, op="="),
    ]

    # Invalid string raises ValueError
    with pytest.raises(ValueError, match="Invalid CIGAR string"):
        parse_cigar("100M_invalid")


def test_cigar_query_and_target_consumption():
    """Verify query and target consumption for M, indels, and diff ops."""
    # Exact match 15643M
    assert cigar_query_consumption("15643M") == 15643
    assert cigar_target_consumption("15643M") == 15643

    # Insertion: 100M2I50M
    # Query consumes M + I = 100 + 2 + 50 = 152
    # Target consumes M = 100 + 50 = 150
    assert cigar_query_consumption("100M2I50M") == 152
    assert cigar_target_consumption("100M2I50M") == 150

    # Deletion: 200M3D10M
    # Query consumes M = 200 + 10 = 210
    # Target consumes M + D = 200 + 3 + 10 = 213
    assert cigar_query_consumption("200M3D10M") == 210
    assert cigar_target_consumption("200M3D10M") == 213

    # Mixed with =, X, I, D
    ops = parse_cigar("10=2X3I4D")
    # Query: 10 + 2 + 3 = 15
    # Target: 10 + 2 + 4 = 16
    assert cigar_query_consumption(ops) == 15
    assert cigar_target_consumption(ops) == 16


def test_cigar_overlap_length():
    """Verify robust calculation of overlap length."""
    assert cigar_overlap_length("15643M") == 15643
    assert cigar_overlap_length("100M2I50M") == 152
    assert cigar_overlap_length("200M3D10M") == 213


def test_calculate_excision_offsets_exact_match():
    """Verify excision offsets for exact 15643M overlap."""
    from_len = 20000
    to_len = 25000
    cigar = "15643M"

    from_off, to_off = calculate_excision_offsets(cigar, from_len, to_len)
    assert from_off == 20000 - 15643  # 4357
    assert to_off == 15643


def test_calculate_excision_offsets_indels():
    """Verify excision offsets for indels (100M2I50M and 200M3D10M)."""
    # 100M2I50M: target consumption 150, query consumption 152
    from_len = 1000
    to_len = 1200
    from_off, to_off = calculate_excision_offsets("100M2I50M", from_len, to_len)
    assert from_off == 1000 - 150  # 850
    assert to_off == 152

    # 200M3D10M: target consumption 213, query consumption 210
    from_off_d, to_off_d = calculate_excision_offsets("200M3D10M", from_len, to_len)
    assert from_off_d == 1000 - 213  # 787
    assert to_off_d == 210


def test_calculate_excision_offsets_validation():
    """Verify input validation: negative lengths or consumption exceeding length."""
    # Negative lengths
    with pytest.raises(ValueError, match="Segment lengths must be non-negative"):
        calculate_excision_offsets("100M", -1, 500)

    with pytest.raises(ValueError, match="Segment lengths must be non-negative"):
        calculate_excision_offsets("100M", 500, -1)

    # Target consumption exceeds from_len
    with pytest.raises(ValueError, match="Target consumption .* exceeds from_len"):
        calculate_excision_offsets("100M", 50, 500)

    # Query consumption exceeds to_len
    with pytest.raises(ValueError, match="Query consumption .* exceeds to_len"):
        calculate_excision_offsets("100M", 500, 50)
