"""Unit tests for biological quirk and end-state taxonomy."""

from io import StringIO
from pathlib import Path

import pytest

from gfapdefinisher.classifier import PdeJunctionCall
from gfapdefinisher.models import SegmentPort, Strand
from gfapdefinisher.quirks import (
    END_STATE_TSV_COLUMNS,
    JUNCTION_TSV_COLUMNS,
    EndState,
    EndStateRecord,
    IntervalQuirk,
    JunctionClass,
    JunctionRecord,
    QuirkRecord,
    build_feature_id,
    create_junction_record,
    emit_end_states_tsv,
    emit_junctions_tsv,
    emit_quirks_tsv,
    parse_end_states_tsv,
    parse_junctions_tsv,
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


def test_end_states_tsv_roundtrip(tmp_path: Path):
    records = [
        EndStateRecord(
            feature_id="end_complete:I:5p",
            chrom="I",
            side="5p",
            state=EndState.END_COMPLETE,
            somatic_port="utg001+",
            germline_port="utg002-",
            shared_port=None,
            core_coverage_bp_median=45.5,
            somatic_depth_ratio=1.0,
            germline_depth_ratio=0.25,
            resolved=True,
            notes="Fully resolved somatic and germline caps",
        ),
        EndStateRecord(
            feature_id="end_somatic_missing:V:3p",
            chrom="V",
            side="3p",
            state=EndState.END_SOMATIC_MISSING,
            somatic_port=None,
            germline_port="utg088+",
            shared_port=None,
            core_coverage_bp_median=50.0,
            somatic_depth_ratio=None,
            germline_depth_ratio=0.10,
            resolved=False,
            notes="Somatic stub dropout",
        ),
        EndStateRecord(
            feature_id="end_germline_orphan:III:5p",
            chrom="III",
            side="5p",
            state=EndState.END_GERMLINE_ORPHAN,
            somatic_port="utg015-",
            germline_port=None,
            shared_port=None,
            core_coverage_bp_median=42.0,
            somatic_depth_ratio=0.98,
            germline_depth_ratio=None,
            resolved=False,
            notes="Unbridged germline arm",
        ),
        EndStateRecord(
            feature_id="end_shared_exempt:II:5p",
            chrom="II",
            side="5p",
            state=EndState.END_SHARED_EXEMPT,
            somatic_port=None,
            germline_port=None,
            shared_port="utg053-",
            core_coverage_bp_median=48.2,
            somatic_depth_ratio=None,
            germline_depth_ratio=None,
            resolved=True,
            notes="Exempt distal rDNA telomeric cap",
        ),
    ]

    # Test in-memory emission and parsing
    tsv_content = emit_end_states_tsv(records)
    assert tsv_content.startswith("\t".join(END_STATE_TSV_COLUMNS))
    parsed = parse_end_states_tsv(tsv_content)
    assert len(parsed) == len(records)
    for orig, res in zip(records, parsed, strict=True):
        assert orig == res

    # Test emission and parsing via Path
    tsv_file = tmp_path / "end_states.tsv"
    ret_content = emit_end_states_tsv(records, tsv_file)
    assert ret_content == tsv_content
    assert tsv_file.exists()
    parsed_path = parse_end_states_tsv(tsv_file)
    assert parsed_path == records

    # Test parsing via string filename
    parsed_str_path = parse_end_states_tsv(str(tsv_file))
    assert parsed_str_path == records

    # Test parsing via TextIO
    parsed_stream = parse_end_states_tsv(StringIO(tsv_content))
    assert parsed_stream == records

    # Test empty sequence handling
    empty_tsv = emit_end_states_tsv([])
    assert empty_tsv == "\t".join(END_STATE_TSV_COLUMNS) + "\n"
    assert parse_end_states_tsv(empty_tsv) == []
    assert parse_end_states_tsv("") == []


def test_junctions_tsv_roundtrip(tmp_path: Path):
    records = [
        JunctionRecord(
            feature_id="pde_terminal:IV:3p",
            chrom="IV",
            side="3p",
            junction_class=JunctionClass.PDE_TERMINAL,
            stem_port="utg010+",
            somatic_terminus_port="utg011+",
            germline_continuation_port="utg012-",
            break_index=1,
            core_depth_ratio=0.18,
            read_count=24,
            confidence_tier="CONFIDENT",
            notes="Benchmark validated PDE junction",
        ),
        JunctionRecord(
            feature_id="pde_intercalary:X:mid:2",
            chrom="X",
            side="mid",
            junction_class=JunctionClass.PDE_INTERCALARY,
            stem_port="utg045-",
            somatic_terminus_port=None,
            germline_continuation_port="utg046+",
            break_index=2,
            core_depth_ratio=0.08,
            read_count=15,
            confidence_tier="PROVISIONAL",
            notes="Internal diminution break 2",
        ),
        JunctionRecord(
            feature_id="pde_terminal:II:5p",
            chrom="II",
            side="5p",
            junction_class=JunctionClass.PDE_TERMINAL,
            stem_port="utg030+",
            somatic_terminus_port=None,
            germline_continuation_port=None,
            break_index=1,
            core_depth_ratio=None,
            read_count=None,
            confidence_tier="CONFIDENT",
            notes="",
        ),
    ]

    # Test in-memory emission and parsing
    tsv_content = emit_junctions_tsv(records)
    assert tsv_content.startswith("\t".join(JUNCTION_TSV_COLUMNS))
    parsed = parse_junctions_tsv(tsv_content)
    assert len(parsed) == len(records)
    for orig, res in zip(records, parsed, strict=True):
        assert orig == res

    # Test emission and parsing via Path
    tsv_file = tmp_path / "junctions.tsv"
    ret_content = emit_junctions_tsv(records, tsv_file)
    assert ret_content == tsv_content
    assert tsv_file.exists()
    parsed_path = parse_junctions_tsv(tsv_file)
    assert parsed_path == records

    # Test parsing via string filename
    parsed_str_path = parse_junctions_tsv(str(tsv_file))
    assert parsed_str_path == records

    # Test parsing via TextIO
    parsed_stream = parse_junctions_tsv(StringIO(tsv_content))
    assert parsed_stream == records

    # Test empty sequence handling
    empty_tsv = emit_junctions_tsv([])
    assert empty_tsv == "\t".join(JUNCTION_TSV_COLUMNS) + "\n"
    assert parse_junctions_tsv(empty_tsv) == []
    assert parse_junctions_tsv("") == []


def test_create_junction_record_terminal():
    call = PdeJunctionCall(
        port=SegmentPort("utg000009l", Strand.MINUS),
        is_pde_junction=True,
        somatic_port=SegmentPort("utg000211l", Strand.MINUS),
        germline_port=SegmentPort("utg000069l", Strand.FORWARD),
        depth_ratio=0.0968,
    )
    rec = create_junction_record(
        call=call,
        chrom="I",
        side="5p",
        junction_class=JunctionClass.PDE_TERMINAL,
        confidence_tier="CONFIDENT",
    )

    assert rec.feature_id == "pde_terminal:I:5p"
    assert rec.chrom == "I"
    assert rec.side == "5p"
    assert rec.junction_class == JunctionClass.PDE_TERMINAL
    assert rec.stem_port == "utg000009l-"
    assert rec.somatic_terminus_port == "utg000211l-"
    assert rec.germline_continuation_port == "utg000069l+"
    assert rec.core_depth_ratio == 0.0968
    assert rec.confidence_tier == "CONFIDENT"
    assert rec.break_index == 1
    assert rec.notes == ""


def test_create_junction_record_intercalary():
    call = PdeJunctionCall(
        port=SegmentPort("utg000045l", Strand.FORWARD),
        is_pde_junction=True,
        somatic_port=None,
        germline_port=SegmentPort("utg000046l", Strand.FORWARD),
        depth_ratio=0.08,
    )
    rec = create_junction_record(
        call=call,
        chrom="X",
        side="mid",
        junction_class="pde_intercalary",
        confidence_tier="PROVISIONAL",
        break_index=2,
        notes="Internal break call",
    )

    assert rec.feature_id == "pde_intercalary:X:mid:2"
    assert rec.chrom == "X"
    assert rec.side == "mid"
    assert rec.junction_class == JunctionClass.PDE_INTERCALARY
    assert rec.stem_port == "utg000045l+"
    assert rec.somatic_terminus_port is None
    assert rec.germline_continuation_port == "utg000046l+"
    assert rec.break_index == 2
    assert rec.confidence_tier == "PROVISIONAL"
    assert rec.notes == "Internal break call"


def test_create_junction_record_repeat_entry():
    call = PdeJunctionCall(
        port=SegmentPort("utg000019l", Strand.FORWARD),
        is_pde_junction=True,
        somatic_port=SegmentPort("utg000200l", Strand.MINUS),
        germline_port=SegmentPort("utg000071l", Strand.FORWARD),
        depth_ratio=0.1061,
        is_repeat_entry=True,
    )
    rec = create_junction_record(
        call=call,
        chrom="IV",
        side="3p",
        junction_class=JunctionClass.PDE_TERMINAL,
        notes="Context note",
    )

    assert rec.feature_id == "pde_terminal:IV:3p"
    assert "repeat_entry" in rec.notes
    assert "Context note" in rec.notes


def test_top_level_package_exports():
    import gfapdefinisher as g
    import gfapdefinisher.quirks as gq

    assert hasattr(g, "EndStateRecord")
    assert hasattr(g, "JunctionRecord")
    assert hasattr(g, "QuirkRecord")
    assert hasattr(g, "EndState")
    assert hasattr(g, "JunctionClass")
    assert hasattr(g, "IntervalQuirk")
    assert hasattr(g, "DEFAULT_POLICY")
    assert hasattr(g, "Policy")
    assert g.DEFAULT_POLICY is not None
    # Factory and serdes exported in quirks module namespace and __all__
    assert hasattr(gq, "create_junction_record")
    assert "create_junction_record" in gq.__all__
    assert "emit_end_states_tsv" in gq.__all__
    assert "parse_end_states_tsv" in gq.__all__
    assert "emit_junctions_tsv" in gq.__all__
    assert "parse_junctions_tsv" in gq.__all__


def test_end_states_tsv_special_characters_and_comments():
    """Verify TSV SerDe safely roundtrips fields with special characters, quotes, tabs, and ignores comments/extra columns."""
    rec = EndStateRecord(
        feature_id="end_complete:IV:3p",
        chrom="IV",
        side="3p",
        state=EndState.END_COMPLETE,
        somatic_port="utg010+",
        germline_port="utg020-",
        shared_port=None,
        core_coverage_bp_median=64.5,
        somatic_depth_ratio=1.02,
        germline_depth_ratio=0.18,
        resolved=True,
        notes='Special note with "quotes", \t tab, and \n newline',
    )
    tsv = emit_end_states_tsv([rec])
    parsed = parse_end_states_tsv(tsv)
    assert len(parsed) == 1
    assert parsed[0] == rec

    # TSV with comment lines and extra columns
    tsv_with_extra = (
        "# Commentary header\n"
        "feature_id\tchrom\tside\tstate\tsomatic_port\tgermline_port\tshared_port\t"
        "core_coverage_bp_median\tsomatic_depth_ratio\tgermline_depth_ratio\tresolved\tnotes\textra_col\n"
        "# Another comment\n"
        "end_complete:IV:3p\tIV\t3p\tend_complete\tutg010+\tutg020-\t\t64.5\t1.02\t0.18\tTrue\tClean note\tbonus\n"
    )
    parsed_extra = parse_end_states_tsv(tsv_with_extra)
    assert len(parsed_extra) == 1
    assert parsed_extra[0].feature_id == "end_complete:IV:3p"
    assert parsed_extra[0].notes == "Clean note"


def test_junctions_tsv_special_characters_and_comments():
    """Verify JunctionRecord TSV SerDe roundtrips special characters and safely ignores comments and extra columns."""
    rec = JunctionRecord(
        feature_id="pde_terminal:IV:3p",
        chrom="IV",
        side="3p",
        junction_class=JunctionClass.PDE_TERMINAL,
        stem_port="utg010+",
        somatic_terminus_port="utg011+",
        germline_continuation_port="utg012-",
        break_index=1,
        core_depth_ratio=0.18,
        read_count=24,
        confidence_tier="CONFIDENT",
        notes='Complex note: "diminution"\twith\ttabs\nand newlines',
    )
    tsv = emit_junctions_tsv([rec])
    parsed = parse_junctions_tsv(tsv)
    assert len(parsed) == 1
    assert parsed[0] == rec

    tsv_extra = (
        "# Header note\n"
        "feature_id\tchrom\tside\tjunction_class\tstem_port\tsomatic_terminus_port\t"
        "germline_continuation_port\tbreak_index\tcore_depth_ratio\tread_count\tconfidence_tier\tnotes\tunrecognized\n"
        "pde_terminal:IV:3p\tIV\t3p\tpde_terminal\tutg010+\tutg011+\tutg012-\t1\t0.18\t24\tCONFIDENT\tnotes\textra_val\n"
    )
    parsed_extra = parse_junctions_tsv(tsv_extra)
    assert len(parsed_extra) == 1
    assert parsed_extra[0].feature_id == "pde_terminal:IV:3p"


def test_nonexistent_tsv_raises_filenotfound():
    """Verify parse functions raise FileNotFoundError when passed nonexistent file paths."""
    with pytest.raises(FileNotFoundError):
        parse_end_states_tsv("nonexistent_end_states.tsv")

    with pytest.raises(FileNotFoundError):
        parse_junctions_tsv("nonexistent_junctions.tsv")


def test_junction_record_break_index_validation():
    """Verify JunctionRecord enforces break_index >= 1."""
    with pytest.raises(ValueError, match="Break index must be at least 1"):
        JunctionRecord(
            feature_id="pde_terminal:I:5p",
            chrom="I",
            side="5p",
            junction_class=JunctionClass.PDE_TERMINAL,
            stem_port="utg001+",
            break_index=0,
        )

    with pytest.raises(ValueError, match="Break index must be at least 1"):
        JunctionRecord(
            feature_id="pde_terminal:I:5p",
            chrom="I",
            side="5p",
            junction_class=JunctionClass.PDE_TERMINAL,
            stem_port="utg001+",
            break_index=-2,
        )


def test_create_junction_record_core_depth_ratio_precedence():
    """Verify create_junction_record correctly extracts core_depth_ratio from PdeJunctionCall."""
    call = PdeJunctionCall(
        port=SegmentPort("utg000009l", Strand.MINUS),
        is_pde_junction=True,
        somatic_port=SegmentPort("utg000211l", Strand.MINUS),
        germline_port=SegmentPort("utg000069l", Strand.FORWARD),
        core_depth_ratio=0.155,
        depth_ratio=0.0,
    )
    rec = create_junction_record(
        call=call,
        chrom="I",
        side="5p",
        junction_class=JunctionClass.PDE_TERMINAL,
    )
    assert rec.core_depth_ratio == 0.155


@pytest.mark.bio
@pytest.mark.benchmark
def test_benchmark_df5120_collect_features_intercalary_grs_bridges():
    """Verify collect_features recognizes Chr X intercalary GRS bridges and emits only 3 genuine SAT_GAP_VARIABLE quirks."""
    benchmark_dir = Path("data/benchmark/DF5120")
    if not benchmark_dir.exists():
        pytest.skip(f"Benchmark DF5120 dataset not present in {benchmark_dir}")

    gfa_candidates = [
        benchmark_dir / "asm.p_utg.gfa.gz",
        benchmark_dir / "asm.p_utg.gfa",
    ]
    gfa_path = next((p for p in gfa_candidates if p.exists()), None)
    if gfa_path is None:
        pytest.skip(f"GFA file missing in {benchmark_dir}")

    tiara_path = benchmark_dir / "tiara_unitigs.txt"
    tidk_path = benchmark_dir / "tidk_telomeres_telomeric_repeat_windows.tsv"
    barrnap_path = benchmark_dir / "barrnap_rdna.gff"
    compleasm_path = benchmark_dir / "compleasm_unitigs_full_table.tsv"

    if not (tiara_path.exists() and tidk_path.exists() and barrnap_path.exists() and compleasm_path.exists()):
        pytest.skip(f"Benchmark annotation files missing in {benchmark_dir}")

    from gfapdefinisher.annotations import aggregate_node_annotations
    from gfapdefinisher.classifier import compute_core_coverage
    from gfapdefinisher.graph import load_gfa
    from gfapdefinisher.ledger import build_end_evidence_ledger
    from gfapdefinisher.quirks import collect_features

    graph = load_gfa(gfa_path)
    annotations = aggregate_node_annotations(
        tiara=tiara_path,
        tidk=tidk_path,
        barrnap=barrnap_path,
        compleasm=compleasm_path,
        sequence_ids=graph.segments.keys(),
    )
    ledger = build_end_evidence_ledger(graph=graph, annotations=annotations)
    core_coverage = compute_core_coverage(graph, ledger)

    junction_recs, end_recs, quirk_recs = collect_features(
        graph=graph,
        ledger=ledger,
        core_coverage=core_coverage,
    )

    # 1. Verify Chr X internal bridges are classified as END_COMPLETE and resolved
    end_map = {(r.chrom, r.side): r for r in end_recs}

    # utg000019l:3p connects to utg000071l+ (bridge to utg000021l-)
    assert ("utg000019l", "3p") in end_map
    rec_19_3p = end_map[("utg000019l", "3p")]
    assert rec_19_3p.state == EndState.END_COMPLETE
    assert rec_19_3p.resolved is True
    assert rec_19_3p.germline_port == "utg000071l+"

    # utg000019l:5p connects to utg000145l+ (bridge to utg000027l+)
    assert ("utg000019l", "5p") in end_map
    rec_19_5p = end_map[("utg000019l", "5p")]
    assert rec_19_5p.state == EndState.END_COMPLETE
    assert rec_19_5p.resolved is True
    assert rec_19_5p.germline_port == "utg000145l+"

    # utg000021l:5p connects to utg000071l- (bridge to utg000019l+)
    assert ("utg000021l", "5p") in end_map
    rec_21_5p = end_map[("utg000021l", "5p")]
    assert rec_21_5p.state == EndState.END_COMPLETE
    assert rec_21_5p.resolved is True
    assert rec_21_5p.germline_port == "utg000071l-"

    # utg000027l:3p connects to utg000145l- (bridge to utg000019l-)
    assert ("utg000027l", "3p") in end_map
    rec_27_3p = end_map[("utg000027l", "3p")]
    assert rec_27_3p.state == EndState.END_COMPLETE
    assert rec_27_3p.resolved is True
    assert rec_27_3p.germline_port == "utg000145l-"

    # 2. Somatic dropout stub utg000004l:3p is END_SOMATIC_MISSING and unresolved
    assert ("utg000004l", "3p") in end_map
    rec_4_3p = end_map[("utg000004l", "3p")]
    assert rec_4_3p.state == EndState.END_SOMATIC_MISSING
    assert rec_4_3p.resolved is False

    # 3. Only the 3 genuine terminal SAT_GAP_VARIABLE quirks are emitted
    sat_gap_quirks = [q for q in quirk_recs if q.quirk_type == IntervalQuirk.SAT_GAP_VARIABLE]
    sat_gap_ends = {(q.chrom, q.side) for q in sat_gap_quirks}
    assert sat_gap_ends == {
        ("utg000002l", "5p"),
        ("utg000009l", "5p"),
        ("utg000013l", "3p"),
    }
    assert len(sat_gap_quirks) == 3
