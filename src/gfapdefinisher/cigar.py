"""CIGAR alignment parser and dovetail excision coordinate engine.

Provides lossless parsing for all standard CIGAR operations (M, =, X, I, D, S, etc.)
and computes exact sequence consumptions and excision offsets for dovetail merges.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

_CIGAR_OP_REGEX = re.compile(r"(\d+)([MIDNSHPX=])", re.IGNORECASE)

QUERY_CONSUMING_OPS = frozenset({"M", "=", "X", "I"})
TARGET_CONSUMING_OPS = frozenset({"M", "=", "X", "D"})


@dataclass(slots=True, frozen=True)
class CigarOp:
    """An individual CIGAR alignment operation."""

    length: int
    op: str


def parse_cigar(cigar_str: str) -> list[CigarOp]:
    """Parse a CIGAR string into a list of CigarOp records.

    Handles all CIGAR operations including match/mismatch (M, =, X),
    indels (I, D), and clipping/skips (S, H, N, P).

    Args:
        cigar_str: CIGAR representation (e.g. '15643M', '100M2I50M', '*').

    Returns:
        List of CigarOp instances in alignment order.

    Raises:
        ValueError: If the string contains unrecognized or malformed tokens.
    """
    if not cigar_str or cigar_str == "*":
        return []

    ops: list[CigarOp] = []
    total_matched = 0
    for match in _CIGAR_OP_REGEX.finditer(cigar_str):
        total_matched += len(match.group(0))
        op_len = int(match.group(1))
        op_char = match.group(2).upper()
        ops.append(CigarOp(length=op_len, op=op_char))

    if total_matched != len(cigar_str):
        raise ValueError(f"Invalid CIGAR string contains unrecognized tokens: {cigar_str!r}")

    return ops


def format_cigar(ops: list[CigarOp]) -> str:
    """Format a list of CigarOp items into a standard CIGAR string."""
    return "".join(f"{op.length}{op.op}" for op in ops)


def cigar_query_consumption(ops: list[CigarOp] | str) -> int:
    """Compute total bases consumed on the query sequence by the CIGAR.

    Query-consuming operations in dovetail overlaps: M, =, X, I.
    """
    parsed = parse_cigar(ops) if isinstance(ops, str) else ops
    return sum(op.length for op in parsed if op.op in QUERY_CONSUMING_OPS)


def cigar_target_consumption(ops: list[CigarOp] | str) -> int:
    """Compute total bases consumed on the target sequence by the CIGAR.

    Target-consuming operations in dovetail overlaps: M, =, X, D.
    """
    parsed = parse_cigar(ops) if isinstance(ops, str) else ops
    return sum(op.length for op in parsed if op.op in TARGET_CONSUMING_OPS)


def cigar_overlap_length(cigar_str: str | list[CigarOp]) -> int:
    """Robustly compute the effective overlap length represented by a CIGAR string.

    For exact kM matches, returns k. For alignments containing indels, returns
    the maximum span consumed between target and query sequences.
    """
    ops = parse_cigar(cigar_str) if isinstance(cigar_str, str) else cigar_str
    t = cigar_target_consumption(ops)
    q = cigar_query_consumption(ops)
    return max(t, q)


def calculate_excision_offsets(
    cigar_str: str | list[CigarOp],
    from_len: int,
    to_len: int,
) -> tuple[int, int]:
    """Calculate exact excision offsets for a dovetail overlap merge.

    In GFA dovetail semantics between from_seg and to_seg:
    - from_seg is the target sequence, overlapping at its 3' end (suffix).
    - to_seg is the query sequence, overlapping at its 5' end (prefix).
    - Target consumption (T) is excised from the 3' end of from_seg.
    - Query consumption (Q) is excised from the 5' end of to_seg.

    Returns:
        (from_offset, to_offset):
        - from_offset = from_len - T (0-based end coordinate of non-overlapping from_seg).
        - to_offset = Q (0-based start coordinate of non-overlapping to_seg).

    Raises:
        ValueError: If segment lengths are negative or overlap consumption exceeds segment length.
    """
    ops = parse_cigar(cigar_str) if isinstance(cigar_str, str) else cigar_str
    t_len = cigar_target_consumption(ops)
    q_len = cigar_query_consumption(ops)

    if from_len < 0 or to_len < 0:
        raise ValueError(f"Segment lengths must be non-negative: from_len={from_len}, to_len={to_len}")
    if t_len > from_len:
        raise ValueError(f"Target consumption ({t_len}) exceeds from_len ({from_len})")
    if q_len > to_len:
        raise ValueError(f"Query consumption ({q_len}) exceeds to_len ({to_len})")

    return (from_len - t_len, q_len)
