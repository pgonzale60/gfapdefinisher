"""Tests for policy parameter governance and anti-magic-number guardrails."""

import ast
from pathlib import Path

import pytest

from gfapdefinisher.policy import DEFAULT_POLICY, PolicyValue

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
    assert DEFAULT_POLICY.stoichiometry.min_germline_depth_ratio.value == 0.04
    assert (
        DEFAULT_POLICY.stoichiometry.min_germline_depth_ratio.value
        < DEFAULT_POLICY.stoichiometry.max_germline_depth_ratio.value
    )


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

    Permits structural 0, 1, -1 and exit codes in cli.py, but flags other numeric constants.
    """
    violations = []
    allowed_numbers = {0, 1, -1, 2}  # Basic indexing, boolean-like offsets, standard exit codes

    for py_file in SRC_DIR.rglob("*.py"):
        if "policy" in py_file.parts:
            continue  # Policy directory is the authorized location for declared values

        tree = ast.parse(py_file.read_text(encoding="utf-8"), filename=str(py_file))
        for node in ast.walk(tree):
            if isinstance(node, ast.Constant) and isinstance(node.value, int | float):
                if node.value not in allowed_numbers and not isinstance(node.value, bool):
                    violations.append(
                        f"{py_file.name}:{node.lineno} has hardcoded number {node.value}. "
                        "Must be read from policy.py instead."
                    )

    assert not violations, "Magic numbers detected in source code:\n" + "\n".join(violations)
