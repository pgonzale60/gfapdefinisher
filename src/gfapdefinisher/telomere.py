"""Native terminal telomere tract scanner for assembly graph sequences.

Provides exact sequence-level scanning for canonical telomeric repeat tracts at
segment physical termini (5' prefix and 3' suffix), computing exact distance/offset
from physical tips, repeat counts, and tract lengths without sliding-window
boundary discretization artifacts.
"""

from __future__ import annotations

from dataclasses import dataclass

from gfapdefinisher.policy import DEFAULT_POLICY


@dataclass(slots=True, frozen=True)
class TerminalTelomereCall:
    """Exact telomeric tract metrics at a sequence physical terminus.

    Attributes:
        tract_length_bp: Span in base pairs from start of first repeat to end of last repeat.
        repeat_count: Number of canonical repeat hexamers chained within the tract.
        offset_bp: Physical distance from the terminus tip to the start of the repeat tract.
        is_terminal_cap: True if offset <= max_terminal_offset_bp and repeat_count >= min_repeat_count.
        motif: Canonical motif queried.
    """

    tract_length_bp: int
    repeat_count: int
    offset_bp: int
    is_terminal_cap: bool
    motif: str = ""


TerminalTelomereResult = TerminalTelomereCall


def scan_terminal_tract(
    sequence: str,
    motif: str,
    terminal_scan_bp: int = DEFAULT_POLICY.telomere.terminal_scan_bp.value,
    max_terminal_offset_bp: int = DEFAULT_POLICY.telomere.max_terminal_offset_bp.value,
    max_internal_gap_bp: int = DEFAULT_POLICY.telomere.max_internal_gap_bp.value,
    min_repeat_count: int = DEFAULT_POLICY.telomere.min_terminal_repeat_count.value,
    is_suffix: bool = False,
) -> TerminalTelomereCall:
    """Scan a sequence terminus for a chained telomeric repeat tract.

    Args:
        sequence: Nucleotide sequence to scan.
        motif: Canonical telomeric repeat motif (e.g. TTAGGC for 3' or GCCTAA for 5').
        terminal_scan_bp: Maximum sequence length inspected from the physical terminus.
        max_terminal_offset_bp: Maximum physical distance from tip to tract for bona fide cap.
        max_internal_gap_bp: Maximum non-repeat gap allowed when chaining adjacent repeats.
        min_repeat_count: Minimum repeat count required to qualify as an authentic cap.
        is_suffix: If True, inspects 3' sequence suffix; if False, inspects 5' prefix.

    Returns:
        TerminalTelomereCall with tract span, repeat count, offset from tip, and cap status.
    """
    if isinstance(terminal_scan_bp, bool):
        is_suffix = terminal_scan_bp
        terminal_scan_bp = DEFAULT_POLICY.telomere.terminal_scan_bp.value

    if not sequence or not motif or len(sequence) < len(motif):
        return TerminalTelomereCall(
            tract_length_bp=0,
            repeat_count=0,
            offset_bp=0,
            is_terminal_cap=False,
            motif=motif,
        )

    seq_len = len(sequence)
    if is_suffix:
        subseq = sequence[-terminal_scan_bp:] if seq_len > terminal_scan_bp else sequence
    else:
        subseq = sequence[:terminal_scan_bp] if seq_len > terminal_scan_bp else sequence

    subseq_upper = subseq.upper()
    motif_upper = motif.upper()
    subseq_len = len(subseq_upper)
    motif_len = len(motif_upper)

    matches: list[int] = []
    pos = 0
    max_pos = subseq_len - motif_len
    while pos <= max_pos:
        idx = subseq_upper.find(motif_upper, pos)
        if idx == -1:
            break
        matches.append(idx)
        pos = idx + motif_len

    if not matches:
        return TerminalTelomereCall(
            tract_length_bp=0,
            repeat_count=0,
            offset_bp=0,
            is_terminal_cap=False,
            motif=motif,
        )

    tracts: list[dict[str, int]] = []
    current_start = matches[0]
    current_end = matches[0] + motif_len
    current_count = 1

    for idx in matches[1:]:
        gap = idx - current_end
        if gap <= max_internal_gap_bp:
            current_end = idx + motif_len
            current_count += 1
        else:
            offset = (subseq_len - current_end) if is_suffix else current_start
            tracts.append(
                {
                    "start": current_start,
                    "end": current_end,
                    "count": current_count,
                    "length": current_end - current_start,
                    "offset": offset,
                }
            )
            current_start = idx
            current_end = idx + motif_len
            current_count = 1

    offset = (subseq_len - current_end) if is_suffix else current_start
    tracts.append(
        {
            "start": current_start,
            "end": current_end,
            "count": current_count,
            "length": current_end - current_start,
            "offset": offset,
        }
    )

    def _tract_priority(t: dict[str, int]) -> tuple[int, int, int]:
        is_cap = (t["offset"] <= max_terminal_offset_bp) and (t["count"] >= min_repeat_count)
        has_min = t["count"] >= min_repeat_count
        category = 2 if is_cap else (1 if has_min else 0)
        if category > 0:
            return (category, -t["offset"], t["count"])
        return (category, t["count"], -t["offset"])

    best = max(tracts, key=_tract_priority)
    is_cap = (best["offset"] <= max_terminal_offset_bp) and (best["count"] >= min_repeat_count)

    return TerminalTelomereCall(
        tract_length_bp=best["length"],
        repeat_count=best["count"],
        offset_bp=best["offset"],
        is_terminal_cap=is_cap,
        motif=motif,
    )
