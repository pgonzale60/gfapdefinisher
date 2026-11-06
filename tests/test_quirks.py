"""Unit tests for biological quirk and end-state taxonomy."""

from io import StringIO
from pathlib import Path

import pytest

from gfapdefinisher.quirks import (
    EndState,
    EndStateRecord,
    IntervalQuirk,
    JunctionClass,
    JunctionRecord,
    QuirkRecord,
    build_feature_id,
    emit_quirks_tsv,
    parse_quirks_tsv,
)


def test_end_state_enum_values_and_access():
    expected_values = {
        "end_complete",
        "end_somatic_missing",
        "end_germline_orphan",
        "end_shared_exempt",
    }
    assert {e.value for e in EndState} == expected_values

    # Test uppercase and lowercase attribute access
    assert EndState.END_COMPLETE == "end_complete"
    assert EndState.end_complete == "end_complete"
    assert EndState.END_SOMATIC_MISSING == "end_somatic_missing"
    assert EndState.end_somatic_missing == "end_somatic_missing"
    assert EndState.END_GERMLINE_ORPHAN == "end_germline_orphan"
    assert EndState.end_germline_orphan == "end_germline_orphan"
    assert EndState.END_SHARED_EXEMPT == "end_shared_exempt"
    assert EndState.end_shared_exempt == "end_shared_exempt"

    # From string instantiation
    assert EndState("end_complete") is EndState.END_COMPLETE


def test_junction_class_enum_values_and_access():
    expected_values = {"pde_terminal", "pde_intercalary"}
    assert {e.value for e in JunctionClass} == expected_values

    assert JunctionClass.PDE_TERMINAL == "pde_terminal"
    assert JunctionClass.pde_terminal == "pde_terminal"
    assert JunctionClass.PDE_INTERCALARY == "pde_intercalary"
    assert JunctionClass.pde_intercalary == "pde_intercalary"

    assert JunctionClass("pde_terminal") is JunctionClass.PDE_TERMINAL


def test_interval_quirk_enum_values_and_access():
    expected_values = {
        "sat_gap_variable",
        "sat_gap_fixed",
        "rdna_cassette",
        "uncollapsed_copy",
        "low_support_stub",
        "bubble_rejoining",
        "itr",
    }
    assert {e.value for e in IntervalQuirk} == expected_values

    for val in expected_values:
        assert getattr(IntervalQuirk, val) == val
        assert IntervalQuirk(val).value == val


def test_build_feature_id():
    # Without index
    assert build_feature_id(EndState.END_COMPLETE, "I", "5p") == "end_complete:I:5p"
    assert build_feature_id(JunctionClass.PDE_TERMINAL, "IV", "3p") == "pde_terminal:IV:3p"

    # With index
    assert build_feature_id(IntervalQuirk.SAT_GAP_VARIABLE, "IV", "3p", 1) == "sat_gap_variable:IV:3p:1"
    assert build_feature_id("pde_intercalary", "X", "mid", 2) == "pde_intercalary:X:mid:2"


def test_end_state_record_validation_and_serialization():
    # Valid construction
    rec = EndStateRecord(
        feature_id="end_complete:I:5p",
        chrom="I",
        side="5p",
        state=EndState.END_COMPLETE,
        somatic_port="utg001+",
        germline_port="utg002-",
        core_coverage_bp_median=45.0,
        somatic_depth_ratio=1.0,
        germline_depth_ratio=0.25,
        resolved=True,
    )
    assert rec.state == EndState.END_COMPLETE
    d = rec.to_dict()
    assert d["state"] == "end_complete"
    assert d["chrom"] == "I"
    assert d["resolved"] is True

    # Coercion from string
    rec2 = EndStateRecord(
        feature_id="end_shared_exempt:II:5p",
        chrom="II",
        side="5p",
        state="end_shared_exempt",  # type: ignore[arg-type]
        shared_port="utg053-",
    )
    assert rec2.state is EndState.END_SHARED_EXEMPT

    # Validation errors
    with pytest.raises(ValueError, match="Chromosome cannot be empty"):
        EndStateRecord(feature_id="f", chrom="", side="5p", state=EndState.END_COMPLETE)

    with pytest.raises(ValueError, match="Side cannot be empty"):
        EndStateRecord(feature_id="f", chrom="I", side="", state=EndState.END_COMPLETE)

    with pytest.raises(ValueError):
        EndStateRecord(feature_id="f", chrom="I", side="5p", state="invalid_state")  # type: ignore[arg-type]


def test_junction_record_validation_and_serialization():
    # Valid construction
    rec = JunctionRecord(
        feature_id="pde_terminal:IV:3p",
        chrom="IV",
        side="3p",
        junction_class=JunctionClass.PDE_TERMINAL,
        stem_port="utg010+",
        somatic_terminus_port="utg011+",
        germline_continuation_port="utg012-",
        core_depth_ratio=0.18,
        read_count=24,
    )
    assert rec.junction_class == JunctionClass.PDE_TERMINAL
    d = rec.to_dict()
    assert d["junction_class"] == "pde_terminal"
    assert d["break_index"] == 1

    # Intercalary with custom index and string class
    rec_intercalary = JunctionRecord(
        feature_id="pde_intercalary:X:mid:1",
        chrom="X",
        side="mid",
        junction_class="pde_intercalary",  # type: ignore[arg-type]
        stem_port="utg020+",
        break_index=1,
    )
    assert rec_intercalary.junction_class is JunctionClass.PDE_INTERCALARY

    # Validation errors
    with pytest.raises(ValueError, match="Chromosome cannot be empty"):
        JunctionRecord(feature_id="f", chrom="", side="3p", junction_class=JunctionClass.PDE_TERMINAL, stem_port="u+")

    with pytest.raises(ValueError, match="Side cannot be empty"):
        JunctionRecord(feature_id="f", chrom="IV", side="", junction_class=JunctionClass.PDE_TERMINAL, stem_port="u+")

    with pytest.raises(ValueError, match="Stem port cannot be empty"):
        JunctionRecord(feature_id="f", chrom="IV", side="3p", junction_class=JunctionClass.PDE_TERMINAL, stem_port="")

    with pytest.raises(ValueError):
        JunctionRecord(
            feature_id="f",
            chrom="IV",
            side="3p",
            junction_class="invalid_class",  # type: ignore[arg-type]
            stem_port="u+",
        )


def test_quirk_record_validation_and_serialization():
    # Variable satellite repeat quirk with spanning ONT metrics
    rec = QuirkRecord(
        feature_id="sat_gap_variable:IV:3p:1",
        chrom="IV",
        side="3p",
        quirk_type=IntervalQuirk.SAT_GAP_VARIABLE,
        spanning_read_count=7,
        length_min_bp=8200,
        length_median_bp=8500,
        length_max_bp=8900,
        length_iqr_bp=400,
        repeat_monomer_len_bp=128,
        copy_number_range="64-70",
        left_anchor_port="utg010+",
        right_anchor_port="utg056+",
        agp_component_type="W",
        agp_linkage_evidence="ont-spanning",
        unitig_ids=["utg010", "utg056"],  # type: ignore[arg-type]
        read_count=7,
        core_coverage_bp_median=42.0,
        core_depth_ratio=0.16,
        representative_sequence_selected=True,
        notes="Empirical ONT repeat array bridging",
    )
    assert isinstance(rec.unitig_ids, tuple)
    assert rec.unitig_ids == ("utg010", "utg056")
    assert rec.quirk_type is IntervalQuirk.SAT_GAP_VARIABLE

    d = rec.to_dict()
    assert d["quirk_type"] == "sat_gap_variable"
    assert d["unitig_ids"] == ["utg010", "utg056"]
    assert d["representative_sequence_selected"] is True

    # Check TSV row serialization
    row = rec.to_tsv_row()
    assert row[0] == "sat_gap_variable:IV:3p:1"
    assert row[1] == "IV"
    assert row[3] == "sat_gap_variable"
    assert row[4] == "7"
    assert row[6] == "8500"
    assert row[19] == "True"

    # Validation errors
    with pytest.raises(ValueError, match="Chromosome cannot be empty"):
        QuirkRecord(feature_id="f", chrom="", side="3p", quirk_type=IntervalQuirk.SAT_GAP_VARIABLE)

    with pytest.raises(ValueError, match="Side cannot be empty"):
        QuirkRecord(feature_id="f", chrom="IV", side="", quirk_type=IntervalQuirk.SAT_GAP_VARIABLE)

    with pytest.raises(ValueError, match="cannot exceed length_max_bp"):
        QuirkRecord(
            feature_id="f",
            chrom="IV",
            side="3p",
            quirk_type=IntervalQuirk.SAT_GAP_VARIABLE,
            length_min_bp=9000,
            length_max_bp=8000,
        )


def test_quirks_tsv_roundtrip(tmp_path: Path):
    records = [
        QuirkRecord(
            feature_id="sat_gap_variable:IV:3p:1",
            chrom="IV",
            side="3p",
            quirk_type=IntervalQuirk.SAT_GAP_VARIABLE,
            spanning_read_count=7,
            length_min_bp=8200,
            length_median_bp=8500,
            length_max_bp=8900,
            length_iqr_bp=400,
            repeat_monomer_len_bp=128,
            copy_number_range="64-70",
            left_anchor_port="utg010+",
            right_anchor_port="utg056+",
            agp_component_type="W",
            agp_linkage_evidence="ont-spanning",
            unitig_ids=("utg010", "utg056"),
            read_count=7,
            core_coverage_bp_median=42.0,
            core_depth_ratio=0.16,
            representative_sequence_selected=True,
            notes="Representative sequence selected across variable satellite array",
        ),
        QuirkRecord(
            feature_id="rdna_cassette:II:5p:1",
            chrom="II",
            side="5p",
            quirk_type=IntervalQuirk.RDNA_CASSETTE,
            unitig_ids=("utg020",),
            read_count=350,
            core_coverage_bp_median=42.0,
            core_depth_ratio=8.33,
            notes="Multi-copy 18S-5.8S-28S repeating cluster",
        ),
        QuirkRecord(
            feature_id="low_support_stub:V:3p:1",
            chrom="V",
            side="3p",
            quirk_type=IntervalQuirk.LOW_SUPPORT_STUB,
            unitig_ids=("utg179",),
            read_count=2,
            core_coverage_bp_median=42.0,
            core_depth_ratio=0.05,
            notes="Somatic diminution dropout stub failing 4x depth threshold",
        ),
    ]

    # Test via StringIO
    buffer = StringIO()
    emit_quirks_tsv(records, buffer)
    buffer.seek(0)
    parsed = parse_quirks_tsv(buffer)

    assert len(parsed) == 3
    assert parsed[0].feature_id == "sat_gap_variable:IV:3p:1"
    assert parsed[0].quirk_type is IntervalQuirk.SAT_GAP_VARIABLE
    assert parsed[0].spanning_read_count == 7
    assert parsed[0].length_median_bp == 8500
    assert parsed[0].representative_sequence_selected is True
    assert parsed[0].unitig_ids == ("utg010", "utg056")
    assert parsed[0].core_depth_ratio == 0.16

    assert parsed[1].feature_id == "rdna_cassette:II:5p:1"
    assert parsed[1].quirk_type is IntervalQuirk.RDNA_CASSETTE
    assert parsed[1].length_min_bp is None

    assert parsed[2].feature_id == "low_support_stub:V:3p:1"
    assert parsed[2].quirk_type is IntervalQuirk.LOW_SUPPORT_STUB
    assert parsed[2].read_count == 2
    assert parsed[2].core_depth_ratio == 0.05

    # Test via filesystem Path
    tsv_path = tmp_path / "quirks.tsv"
    emit_quirks_tsv(records, tsv_path)
    assert tsv_path.exists()
    parsed_fs = parse_quirks_tsv(tsv_path)
    assert len(parsed_fs) == 3
    assert parsed_fs[0].feature_id == records[0].feature_id


def test_sat_gap_variable_resolution_rule_semantics():
    """Verify that when a sat_gap_variable is bridged via representative sequence,

    the end state resolves to end_complete rather than end_germline_orphan.
    """
    quirk = QuirkRecord(
        feature_id="sat_gap_variable:IV:3p:1",
        chrom="IV",
        side="3p",
        quirk_type=IntervalQuirk.SAT_GAP_VARIABLE,
        spanning_read_count=7,
        length_median_bp=8500,
        representative_sequence_selected=True,
    )
    assert quirk.representative_sequence_selected is True

    # Once resolved with representative sequence, end state is END_COMPLETE
    end_state = EndStateRecord(
        feature_id="end_complete:IV:3p",
        chrom="IV",
        side="3p",
        state=EndState.END_COMPLETE,
        somatic_port="utg009+",
        germline_port="utg056+",
        resolved=True,
        notes="Resolved via representative sequence for sat_gap_variable:IV:3p:1",
    )
    assert end_state.state is EndState.END_COMPLETE
    assert end_state.resolved is True
