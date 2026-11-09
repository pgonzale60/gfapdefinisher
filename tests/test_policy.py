"""Tests for policy parameter governance and anti-magic-number guardrails."""

import ast
from pathlib import Path

import pytest

from gfapdefinisher.policy import DEFAULT_POLICY, PolicyValue
from gfapdefinisher.store import DEFAULT_BATCH_SIZE, MIN_A_LINE_TOKENS

SRC_DIR = Path(__file__).resolve().parent.parent / "src" / "gfapdefinisher"


def test_policy_value_validation():
    # Valid policy value
    pv = PolicyValue(value=42, justification="Documented reason", provenance="Literature")
    assert pv.value == 42
    assert not pv.provisional

    # Empty justification rejected
    with pytest.raises(ValueError, match="must have a non-empty justification"):
        PolicyValue(value=42, justification="", provenance="Literature")

    # Empty provenance rejected
    with pytest.raises(ValueError, match="must have a non-empty provenance"):
        PolicyValue(value=42, justification="Reason", provenance="  ")


def test_default_policy_completeness():
    assert DEFAULT_POLICY.telomere.motif_forward.value == "TTAGGC"
    assert DEFAULT_POLICY.telomere.motif_reverse.value == "GCCTAA"
    assert DEFAULT_POLICY.telomere.min_terminal_repeat_count.value == 30
    assert not hasattr(DEFAULT_POLICY.telomere, "min_terminal_repeat_fraction")
    assert (
        "biological calibration on reference benchmark" in DEFAULT_POLICY.telomere.min_terminal_repeat_count.provenance
    )
    assert DEFAULT_POLICY.telomere.terminal_scan_bp.value == 10000
    assert DEFAULT_POLICY.telomere.terminal_scan_bp.provisional is True
    assert DEFAULT_POLICY.telomere.max_terminal_offset_bp.value == 25
    assert DEFAULT_POLICY.telomere.max_terminal_offset_bp.provisional is True
    assert "biological calibration on reference benchmark" in DEFAULT_POLICY.telomere.max_terminal_offset_bp.provenance
    assert DEFAULT_POLICY.telomere.max_internal_gap_bp.value == 25
    assert DEFAULT_POLICY.telomere.max_internal_gap_bp.provisional is True
    assert "biological calibration on reference benchmark" in DEFAULT_POLICY.telomere.max_internal_gap_bp.provenance
    assert DEFAULT_POLICY.stoichiometry.min_stem_core_ratio.value == 0.50
    assert DEFAULT_POLICY.stoichiometry.min_germline_depth_ratio.value == 0.01
    assert DEFAULT_POLICY.stoichiometry.max_germline_depth_ratio.value == 0.33
    min_germline = DEFAULT_POLICY.stoichiometry.min_germline_depth_ratio.value
    max_germline = DEFAULT_POLICY.stoichiometry.max_germline_depth_ratio.value
    min_stem_core = DEFAULT_POLICY.stoichiometry.min_stem_core_ratio.value
    assert 0.0 < min_germline < max_germline < min_stem_core <= 1.0
    assert DEFAULT_POLICY.pde.max_somatic_cap_hops.value == 2
    assert DEFAULT_POLICY.pde.max_somatic_cap_hops.provisional is False
    assert DEFAULT_POLICY.ont.min_anchor_bp.value == 1000
    assert DEFAULT_POLICY.ont.max_boundary_slack_bp.value == 1500
    assert DEFAULT_POLICY.ont.min_mapq.value == 0
    assert DEFAULT_POLICY.ont.min_spanning_reads.value == 2
    assert DEFAULT_POLICY.ont.quantile_num.value == 4
    assert DEFAULT_POLICY.rdna.max_rdna_unitig_length_bp.value == 200000
    assert DEFAULT_POLICY.rdna.max_rdna_unitig_length_bp.provisional is False


def test_no_species_hardcoding_in_src():
    """Anti-overfitting rule: no source file may hardcode benchmark or held-out species."""
    banned_species = ["DF5120", "JU1382", "PS1017", "PS2068"]
    violations = []

    for py_file in SRC_DIR.rglob("*.py"):
        text = py_file.read_text(encoding="utf-8")
        for sp in banned_species:
            if sp in text:
                violations.append(f"{py_file.relative_to(SRC_DIR)}: hardcoded species '{sp}'")

    assert not violations, "Species hardcoding detected in source code:\n" + "\n".join(violations)


def test_no_magic_numbers_in_src():
    """Anti-heuristic rule: numeric thresholds must be imported from policy, not hardcoded.

    Permits structural 0, 1, -1 and exit codes in cli.py, parsing format column indices and token counts,
    and data store batch sizes, but flags policy-level biological thresholds.
    """
    violations = []
    allowed_numbers = {0, 1, -1, 2}  # Basic indexing, boolean-like offsets, standard exit codes

    # Explicit non-policy operational exemptions:
    # - store.py: DEFAULT_BATCH_SIZE, MIN_A_LINE_TOKENS (parsing token counts & batch buffers)
    # - annotations.py: column indices for standard 0-indexed TSV/GFF/BED format parsing (0..15)
    # - resolver.py: column indices for standard 0-indexed PAF format parsing (0..15)
    exemptions_by_file: dict[str, set[int | float]] = {
        "store.py": {DEFAULT_BATCH_SIZE, MIN_A_LINE_TOKENS},
        "annotations.py": set(range(16)),
        "resolver.py": set(range(16)),
    }

    for py_file in SRC_DIR.rglob("*.py"):
        if "policy" in py_file.parts:
            continue  # Policy directory is the authorized location for declared values

        tree = ast.parse(py_file.read_text(encoding="utf-8"), filename=str(py_file))
        file_exemptions = exemptions_by_file.get(py_file.name, set())
        for node in ast.walk(tree):
            if isinstance(node, ast.Constant) and isinstance(node.value, int | float):
                val = node.value
                if isinstance(val, bool):
                    continue
                if val in allowed_numbers or val in file_exemptions:
                    continue
                violations.append(
                    f"{py_file.name}:{node.lineno} has hardcoded number {val}. Must be read from policy.py instead."
                )

    assert not violations, "Magic numbers detected in source code:\n" + "\n".join(violations)
