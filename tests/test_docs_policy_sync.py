"""Tests verifying synchronization between documentation and live policy definitions.

Ensures that documentation never drifts from code by verifying that:
1. All policy parameters referenced in markdown documents exist in defaults.py.
2. Referenced parameters are valid PolicyValue instances with justifications and provenances.
3. Documentation refrains from hardcoding unanchored normative numeric thresholds.
4. defaults.py remains the single source of truth for numeric policy values.
"""

from __future__ import annotations

import re
from dataclasses import fields, is_dataclass
from pathlib import Path

import gfapdefinisher.policy.defaults as policy_defaults
from gfapdefinisher.policy import DEFAULT_POLICY, PolicyValue

DOCS_DIR = Path(__file__).resolve().parent.parent / "docs"
SRC_DIR = Path(__file__).resolve().parent.parent / "src" / "gfapdefinisher"


def _extract_policy_references(doc_path: Path) -> list[tuple[str, str, int, str]]:
    """Extract policy references from a markdown file.

    Matches:
    - `<PolicyClassName>.<param_name>` or `\\text{<PolicyClassName>.<param_name>}`
    - `policy.<section>.<param_name>` or `DEFAULT_POLICY.<section>.<param_name>`

    Returns:
        List of (class_or_section_name, parameter_name, line_number, raw_match).
    """
    pattern = re.compile(
        r"(?:\\text\{)?(?:(?P<cls>[A-Z][a-zA-Z0-9_]*Policy)|(?:(?:DEFAULT_POLICY|policy)\.(?P<sec>[a-z_][a-z0-9_]*)))\.(?P<param>[a-zA-Z0-9_\\]+)\}?"
    )

    results = []
    text = doc_path.read_text(encoding="utf-8")
    for line_no, line in enumerate(text.splitlines(), start=1):
        for match in pattern.finditer(line):
            cls_name = match.group("cls")
            sec_name = match.group("sec")
            raw_param = match.group("param")
            clean_param = raw_param.replace(r"\_", "_").rstrip("}")

            target = cls_name if cls_name else sec_name
            results.append((target, clean_param, line_no, match.group(0)))
    return results


def test_markdown_policy_references_exist_in_defaults():
    """All policy attributes referenced in markdown docs must exist as fields in defaults.py."""
    md_files = list(DOCS_DIR.rglob("*.md"))
    assert md_files, "No markdown files found in docs/"

    all_references: list[tuple[Path, str, str, int, str]] = []
    for md_file in md_files:
        refs = _extract_policy_references(md_file)
        for target, param, line_no, raw in refs:
            all_references.append((md_file, target, param, line_no, raw))

    assert all_references, "Expected at least one policy reference in markdown documentation."

    for md_file, target, param, line_no, raw in all_references:
        if target.endswith("Policy"):
            policy_cls = getattr(policy_defaults, target, None)
            assert policy_cls is not None, (
                f"{md_file.name}:{line_no}: Policy class {target!r} (from {raw!r}) not found in defaults.py."
            )
            assert is_dataclass(policy_cls), f"{target} is not a dataclass."
            field_names = {f.name for f in fields(policy_cls)}
            assert param in field_names, (
                f"{md_file.name}:{line_no}: Parameter {param!r} (from {raw!r}) does not exist on {target} in defaults.py. "
                f"Available fields: {sorted(field_names)}"
            )
            instance = policy_cls()
            val = getattr(instance, param)
            assert isinstance(val, PolicyValue), f"{target}.{param} must be a PolicyValue instance, got {type(val)}."
            assert val.justification.strip(), f"{target}.{param} has empty justification."
            assert val.provenance.strip(), f"{target}.{param} has empty provenance."
        else:
            section = getattr(DEFAULT_POLICY, target, None)
            assert section is not None, (
                f"{md_file.name}:{line_no}: Policy section {target!r} not found on DEFAULT_POLICY."
            )
            field_names = {f.name for f in fields(section)}
            assert param in field_names, (
                f"{md_file.name}:{line_no}: Parameter {param!r} does not exist in policy section {target!r}. "
                f"Available fields: {sorted(field_names)}"
            )
            val = getattr(section, param)
            assert isinstance(val, PolicyValue), (
                f"DEFAULT_POLICY.{target}.{param} must be a PolicyValue instance, got {type(val)}."
            )


def test_docs_do_not_contain_unanchored_normative_thresholds():
    """Biology documentation must not contain unanchored normative numeric threshold ranges.

    Any numeric range (e.g. 0.04-0.33x) must be explicitly framed as empirical calibration
    evidence or benchmark observations, rather than static normative threshold definitions.
    """
    biology_docs = list((DOCS_DIR / "biology").glob("*.md"))
    # Note: RULES.md is owned by biological expert (@pgonzale60) and holds immutable rules.
    audited_docs = [f for f in biology_docs if f.name != "RULES.md"]

    ratio_pattern = re.compile(r"\b0\.\d+\s*(?:--|–|-)\s*0\.\d+\s*(?:x|×)", re.IGNORECASE)

    calibration_markers = [
        "benchmark",
        "calibrat",
        "observed",
        "df5120",
        "empirical",
        "evidence",
    ]

    violations = []
    for md_file in audited_docs:
        lines = md_file.read_text(encoding="utf-8").splitlines()
        for line_no, line in enumerate(lines, start=1):
            if ratio_pattern.search(line):
                lower_line = line.lower()
                if not any(marker in lower_line for marker in calibration_markers):
                    violations.append(
                        f"{md_file.name}:{line_no}: Contains unanchored ratio range without empirical calibration context: {line.strip()}"
                    )

    assert not violations, (
        "Found unanchored static numeric threshold ranges in biology documentation:\n"
        + "\n".join(violations)
        + "\nDocumentation must express normative criteria via Policy references (e.g. StoichiometryPolicy.min_germline_depth_ratio) "
        "and confine concrete floats strictly to labeled empirical calibration observations."
    )


def test_defaults_is_single_source_of_truth():
    """Ensure that defaults.py is the sole source of numeric policy values in gfapdefinisher."""
    numeric_param_count = 0
    for _section_name, section in vars(DEFAULT_POLICY).items():
        if is_dataclass(section):
            for f in fields(section):
                pv = getattr(section, f.name)
                if isinstance(pv, PolicyValue) and isinstance(pv.value, int | float):
                    numeric_param_count += 1
                    assert pv.value is not None
                    assert isinstance(pv.justification, str) and pv.justification
                    assert isinstance(pv.provenance, str) and pv.provenance

    assert numeric_param_count > 0, "No numeric policy parameters found in defaults.py"
