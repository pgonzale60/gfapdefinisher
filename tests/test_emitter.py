"""Comprehensive unit tests for FASTA and AGP assembly emitter.

Verifies:
- SequenceProvenance enum classifications without subjective confidence tiers.
- PathComponent and AgpRecord dataclasses with empirical risk provenance tags.
- Dovetail overlap-aware sequence stitching (BIO-005).
- FASTA and AGP v2.1 serialization and round-trip parsing.
- QuirkRecord sequence provenance and unpolished flag serialization.
- ONT-derived bridge representative sequence splicing (BIO-004).
"""

from io import StringIO
from pathlib import Path

import pytest

from gfapdefinisher.emitter import (
    AgpGapRow,
    AgpRecord,
    FastaAgpEmitter,
    PathComponent,
    ScaffoldPlan,
    ScaffoldStep,
    SequenceProvenance,
    assemble_path_sequence,
    assemble_scaffold_sequence,
    compute_break_site_summary,
    emit_agp,
    emit_fasta,
    emit_grs_bed,
    emit_scaffold,
    emit_summary_tsv,
    emit_telomere_bed,
    format_summary_table,
    nominate_super_scaffolds,
    parse_agp,
    parse_summary_tsv,
)
from gfapdefinisher.graph import BidirectedGraph
from gfapdefinisher.models import SegmentPort, Strand
from gfapdefinisher.policy import DEFAULT_POLICY
from gfapdefinisher.resolutions import (
    IntervalQuirk,
    JunctionClass,
    JunctionRecord,
    QuirkRecord,
    ResolutionRecord,
    ResolutionType,
    emit_quirks_tsv,
    parse_quirks_tsv,
)
from gfapdefinisher.resolver import OntBridgingEngine


def test_sequence_provenance_enum_values_and_access():
    expected_values = {
        "HiFi-core",
        "HiFi-germline",
        "HiFi-1x",
        "ONT-derived",
        "rDNA-chained",
    }
    assert {p.value for p in SequenceProvenance} == expected_values

    # Uppercase attribute access
    assert SequenceProvenance.HIFI_CORE == "HiFi-core"
    assert SequenceProvenance.HIFI_GERMLINE == "HiFi-germline"
    assert SequenceProvenance.HIFI_1X == "HiFi-1x"
    assert SequenceProvenance.ONT_DERIVED == "ONT-derived"
    assert SequenceProvenance.RDNA_CHAINED == "rDNA-chained"

    # Lowercase attribute access
    assert SequenceProvenance.hifi_core == "HiFi-core"
    assert SequenceProvenance.hifi_germline == "HiFi-germline"
    assert SequenceProvenance.hifi_1x == "HiFi-1x"
    assert SequenceProvenance.ont_derived == "ONT-derived"
    assert SequenceProvenance.rdna_chained == "rDNA-chained"

    # From string instantiation
    assert SequenceProvenance("HiFi-core") is SequenceProvenance.HIFI_CORE
    assert SequenceProvenance("ONT-derived") is SequenceProvenance.ONT_DERIVED


def test_no_confidence_attributes_in_emitter():
    """Verify strictly no 'confidence' field or attribute in emitter classes."""
    for cls in (SequenceProvenance, PathComponent, AgpRecord):
        fields_and_attrs = set(dir(cls))
        if hasattr(cls, "__dataclass_fields__"):
            fields_and_attrs.update(cls.__dataclass_fields__.keys())
        for attr in fields_and_attrs:
            assert "confidence" not in attr.lower(), f"Forbidden 'confidence' found in {cls.__name__}.{attr}"


def test_path_component_validation_and_properties():
    # Valid component
    comp = PathComponent(
        segment_id="utg001",
        strand="+",
        sequence="ACGTACGT",
        overlap_to_next_bp=4,
        sequence_provenance=SequenceProvenance.HIFI_CORE,
        unpolished=False,
    )
    assert comp.segment_id == "utg001"
    assert comp.strand == "+"
    assert comp.overlap_to_next_bp == 4
    assert comp.sequence_provenance is SequenceProvenance.HIFI_CORE
    assert comp.unpolished is False
    assert comp.is_unpolished is False

    # Strand enum coercion
    comp_strand_enum = PathComponent(
        segment_id="utg002",
        strand=Strand.MINUS,
        sequence="TTTTCCCC",
        overlap_to_next_bp=0,
        sequence_provenance="HiFi-germline",
        unpolished=False,
    )
    assert comp_strand_enum.strand == "-"
    assert comp_strand_enum.sequence_provenance is SequenceProvenance.HIFI_GERMLINE

    # String provenance coercion
    comp_str_prov = PathComponent(
        segment_id="read_ont_01",
        strand="-",
        sequence="GATTACA",
        overlap_to_next_bp=0,
        sequence_provenance="ONT-derived",
        unpolished=True,
    )
    assert comp_str_prov.sequence_provenance is SequenceProvenance.ONT_DERIVED
    assert comp_str_prov.unpolished is True
    assert comp_str_prov.is_unpolished is True

    # Validation errors
    with pytest.raises(ValueError, match="segment_id cannot be empty"):
        PathComponent(segment_id="", strand="+", sequence="A")

    with pytest.raises(ValueError, match="Invalid strand"):
        PathComponent(segment_id="u", strand="?", sequence="A")

    with pytest.raises(ValueError, match="overlap_to_next_bp cannot be negative"):
        PathComponent(segment_id="u", strand="+", sequence="A", overlap_to_next_bp=-1)

    with pytest.raises(ValueError, match="cannot exceed sequence length"):
        PathComponent(segment_id="u", strand="+", sequence="ACGT", overlap_to_next_bp=10)


def test_agp_record_validation_and_line_serialization():
    rec = AgpRecord(
        object_id="chrI",
        object_beg=1,
        object_end=1000,
        part_number=1,
        component_type="W",
        component_id="utg001",
        component_beg=1,
        component_end=1000,
        orientation="+",
        sequence_provenance=SequenceProvenance.HIFI_CORE,
        unpolished=False,
    )
    assert rec.is_unpolished is False
    line = rec.to_line()
    assert line == "chrI\t1\t1000\t1\tW\tutg001\t1\t1000\t+\t# provenance=HiFi-core;unpolished=False"

    # Unpolished ONT record with comment
    ont_rec = AgpRecord(
        object_id="chrI",
        object_beg=1001,
        object_end=1500,
        part_number=2,
        component_type="W",
        component_id="ont_read_009",
        component_beg=1,
        component_end=500,
        orientation="+",
        sequence_provenance=SequenceProvenance.ONT_DERIVED,
        unpolished=True,
        comment="BIO-004 satellite gap bridge",
    )
    assert ont_rec.is_unpolished is True
    ont_line = ont_rec.to_line()
    assert (
        ont_line
        == "chrI\t1001\t1500\t2\tW\tont_read_009\t1\t500\t+\t# provenance=ONT-derived;unpolished=True; BIO-004 satellite gap bridge"
    )

    # Gap record (type 'U') without provenance
    gap_rec = AgpRecord(
        object_id="chrI",
        object_beg=1501,
        object_end=1600,
        part_number=3,
        component_type="U",
        component_id="100",
        component_beg="scaffold",
        component_end="yes",
        orientation="unbridged",
        comment="unbridged gap",
    )
    gap_line = gap_rec.to_line()
    assert gap_line == "chrI\t1501\t1600\t3\tU\t100\tscaffold\tyes\tunbridged\t# unbridged gap"

    # Validation errors
    with pytest.raises(ValueError, match="object_id cannot be empty"):
        AgpRecord("", 1, 10, 1, "W", "u", 1, 10, "+")

    with pytest.raises(ValueError, match="object_beg must be >= 1"):
        AgpRecord("chrI", 0, 10, 1, "W", "u", 1, 10, "+")

    with pytest.raises(ValueError, match="cannot be less than object_beg"):
        AgpRecord("chrI", 10, 5, 1, "W", "u", 1, 10, "+")

    with pytest.raises(ValueError, match="part_number must be >= 1"):
        AgpRecord("chrI", 1, 10, 0, "W", "u", 1, 10, "+")


def test_assemble_scaffold_sequence_empty():
    seq, records = assemble_scaffold_sequence([])
    assert seq == ""
    assert records == []


def test_assemble_scaffold_sequence_single_component():
    comp = PathComponent(
        segment_id="utg001",
        strand="+",
        sequence="ACGTACGTACGT",
        overlap_to_next_bp=0,
        sequence_provenance=SequenceProvenance.HIFI_CORE,
    )
    seq, records = assemble_scaffold_sequence([comp], scaffold_id="chrX")
    assert seq == "ACGTACGTACGT"
    assert len(records) == 1
    r = records[0]
    assert r.object_id == "chrX"
    assert r.object_beg == 1
    assert r.object_end == 12
    assert r.part_number == 1
    assert r.component_id == "utg001"
    assert r.component_beg == 1
    assert r.component_end == 12
    assert r.orientation == "+"
    assert r.sequence_provenance is SequenceProvenance.HIFI_CORE
    assert r.unpolished is False


def test_assemble_scaffold_sequence_dovetail_trimming_bio_005():
    """Verify lossless dovetail overlap trimming across oriented components (BIO-005)."""
    # Segment 1 ends with 'TCG' (3 bp overlap)
    # Segment 2 begins with 'TCG', ends with 'AAA' (3 bp overlap)
    # Segment 3 begins with 'AAA'
    c1 = PathComponent(
        segment_id="utg_a",
        strand="+",
        sequence="ATCGATCG",  # 8 bp
        overlap_to_next_bp=3,
        sequence_provenance=SequenceProvenance.HIFI_CORE,
    )
    c2 = PathComponent(
        segment_id="utg_b",
        strand="+",
        sequence="TCGCCCCAAA",  # 10 bp
        overlap_to_next_bp=3,
        sequence_provenance=SequenceProvenance.HIFI_1X,
    )
    c3 = PathComponent(
        segment_id="utg_c",
        strand="+",
        sequence="AAATTTT",  # 7 bp
        overlap_to_next_bp=0,
        sequence_provenance=SequenceProvenance.HIFI_GERMLINE,
    )

    seq, records = assemble_scaffold_sequence([c1, c2, c3], scaffold_id="scaffold_1")

    # Stitched sequence:
    # c1: "ATCGATCG" (8 bp)
    # c2: "CCCCAAA" (7 bp, trimmed 3 bp 'TCG')
    # c3: "TTTT" (4 bp, trimmed 3 bp 'AAA')
    # Total length: 8 + 7 + 4 = 19 bp
    expected_seq = "ATCGATCGCCCCAAATTTT"
    assert seq == expected_seq
    assert len(seq) == 19

    assert len(records) == 3

    r1, r2, r3 = records
    # r1: 1..8
    assert r1.object_beg == 1
    assert r1.object_end == 8
    assert r1.component_beg == 1
    assert r1.component_end == 8
    assert r1.sequence_provenance is SequenceProvenance.HIFI_CORE

    # r2: 9..15 (length 7)
    assert r2.object_beg == 9
    assert r2.object_end == 15
    assert r2.component_beg == 4  # 1-based start after 3 bp trimmed
    assert r2.component_end == 10
    assert r2.sequence_provenance is SequenceProvenance.HIFI_1X

    # r3: 16..19 (length 4)
    assert r3.object_beg == 16
    assert r3.object_end == 19
    assert r3.component_beg == 4  # 1-based start after 3 bp trimmed
    assert r3.component_end == 7
    assert r3.sequence_provenance is SequenceProvenance.HIFI_GERMLINE

    # Verify alias works identically
    alias_seq, alias_records = assemble_path_sequence([c1, c2, c3], scaffold_id="scaffold_1")
    assert alias_seq == seq
    assert len(alias_records) == len(records)


def test_assemble_scaffold_sequence_minus_strand_component():
    """Verify coordinate translation when component is oriented on minus strand."""
    c1 = PathComponent(
        segment_id="utg_pos",
        strand="+",
        sequence="AAAAAGGG",  # 8 bp
        overlap_to_next_bp=3,
        sequence_provenance=SequenceProvenance.HIFI_CORE,
    )
    # Minus strand component: length 10. Trim 3 bp from 5' of oriented sequence.
    # Original component coordinates used should be 1 .. (10 - 3) = 1..7.
    c2 = PathComponent(
        segment_id="utg_neg",
        strand="-",
        sequence="GGGTTTTTTT",  # 10 bp
        overlap_to_next_bp=0,
        sequence_provenance=SequenceProvenance.HIFI_GERMLINE,
    )

    seq, records = assemble_scaffold_sequence([c1, c2], scaffold_id="chrIV")
    expected_seq = "AAAAAGGGTTTTTTT"
    assert seq == expected_seq
    assert len(seq) == 15

    r2 = records[1]
    assert r2.object_beg == 9
    assert r2.object_end == 15
    assert r2.orientation == "-"
    assert r2.component_beg == 1
    assert r2.component_end == 7  # 7 bp used
    assert r2.component_end - r2.component_beg + 1 == 7


def test_assemble_scaffold_sequence_ont_bridge_bio_004():
    """Verify splicing of unpolished ONT-derived bridge representative sequence (BIO-004)."""
    left_hifi = PathComponent(
        segment_id="utg_left",
        strand="+",
        sequence="A" * 100,
        overlap_to_next_bp=0,
        sequence_provenance=SequenceProvenance.HIFI_CORE,
        unpolished=False,
    )
    ont_bridge = PathComponent(
        segment_id="read_ont_999",
        strand="+",
        sequence="C" * 250,  # 250 bp variable satellite array
        overlap_to_next_bp=0,
        sequence_provenance=SequenceProvenance.ONT_DERIVED,
        unpolished=True,
    )
    right_hifi = PathComponent(
        segment_id="utg_right",
        strand="+",
        sequence="T" * 100,
        overlap_to_next_bp=0,
        sequence_provenance=SequenceProvenance.HIFI_CORE,
        unpolished=False,
    )

    seq, records = assemble_scaffold_sequence([left_hifi, ont_bridge, right_hifi], scaffold_id="chrV")
    assert seq == ("A" * 100) + ("C" * 250) + ("T" * 100)
    assert len(seq) == 450

    assert len(records) == 3
    r_bridge = records[1]
    assert r_bridge.object_id == "chrV"
    assert r_bridge.object_beg == 101
    assert r_bridge.object_end == 350
    assert r_bridge.sequence_provenance is SequenceProvenance.ONT_DERIVED
    assert r_bridge.unpolished is True
    assert r_bridge.is_unpolished is True
    assert "# provenance=ONT-derived;unpolished=True" in r_bridge.to_line()


def test_assemble_scaffold_sequence_overlap_exceeds_length_raises():
    c1 = PathComponent(segment_id="u1", strand="+", sequence="AAA", overlap_to_next_bp=3)
    c2 = PathComponent(segment_id="u2", strand="+", sequence="CCC", overlap_to_next_bp=0)
    with pytest.raises(ValueError, match="exceeds or equals component sequence length"):
        assemble_scaffold_sequence([c1, c2])


def test_emit_fasta(tmp_path: Path):
    scaffolds = {
        "chrI": "ATCGATCGATCG",
        "chrII": "GGGGCCCCAAAATTTT",
    }

    # StringIO emission
    buf = StringIO()
    emit_fasta(scaffolds, buf)
    fasta_text = buf.getvalue()
    assert fasta_text == ">chrI\nATCGATCGATCG\n>chrII\nGGGGCCCCAAAATTTT\n"

    # Line wrapping test
    buf_wrap = StringIO()
    emit_fasta({"chrI": "ACGTACGTACGT"}, buf_wrap, line_width=4)
    assert buf_wrap.getvalue() == ">chrI\nACGT\nACGT\nACGT\n"

    # File path emission
    f_path = tmp_path / "assembly.fasta"
    emit_fasta(scaffolds, f_path)
    assert f_path.exists()
    assert f_path.read_text(encoding="utf-8") == fasta_text


def test_emit_and_parse_agp_roundtrip(tmp_path: Path):
    records = [
        AgpRecord(
            object_id="chrI",
            object_beg=1,
            object_end=1000,
            part_number=1,
            component_type="W",
            component_id="utg001",
            component_beg=1,
            component_end=1000,
            orientation="+",
            sequence_provenance=SequenceProvenance.HIFI_CORE,
            unpolished=False,
        ),
        AgpRecord(
            object_id="chrI",
            object_beg=1001,
            object_end=1800,
            part_number=2,
            component_type="W",
            component_id="read_ont_88",
            component_beg=1,
            component_end=800,
            orientation="+",
            sequence_provenance=SequenceProvenance.ONT_DERIVED,
            unpolished=True,
            comment="variable satellite repeat bridge",
        ),
        AgpRecord(
            object_id="chrI",
            object_beg=1801,
            object_end=2800,
            part_number=3,
            component_type="W",
            component_id="utg002",
            component_beg=1,
            component_end=1000,
            orientation="-",
            sequence_provenance=SequenceProvenance.HIFI_GERMLINE,
            unpolished=False,
        ),
        AgpRecord(
            object_id="chrI",
            object_beg=2801,
            object_end=2900,
            part_number=4,
            component_type="U",
            component_id="100",
            component_beg="scaffold",
            component_end="yes",
            orientation="unbridged",
            comment="unbridged gap",
        ),
    ]

    # StringIO round-trip
    buf = StringIO()
    emit_agp(records, buf)
    agp_content = buf.getvalue()
    assert agp_content.startswith("##agp-version 2.1\n")

    parsed = parse_agp(StringIO(agp_content))
    assert len(parsed) == 4

    assert parsed[0].object_id == "chrI"
    assert parsed[0].object_beg == 1
    assert parsed[0].object_end == 1000
    assert parsed[0].sequence_provenance is SequenceProvenance.HIFI_CORE
    assert parsed[0].unpolished is False

    assert parsed[1].object_id == "chrI"
    assert parsed[1].object_beg == 1001
    assert parsed[1].object_end == 1800
    assert parsed[1].sequence_provenance is SequenceProvenance.ONT_DERIVED
    assert parsed[1].unpolished is True
    assert "variable satellite repeat bridge" in (parsed[1].comment or "")

    assert parsed[2].sequence_provenance is SequenceProvenance.HIFI_GERMLINE
    assert parsed[2].orientation == "-"

    assert parsed[3].component_type == "U"
    assert parsed[3].component_beg == "scaffold"
    assert parsed[3].component_end == "yes"

    # File path round-trip
    f_path = tmp_path / "assembly.agp"
    emit_agp(records, f_path)
    parsed_fs = parse_agp(f_path)
    assert len(parsed_fs) == 4
    assert parsed_fs[1].sequence_provenance is SequenceProvenance.ONT_DERIVED


def test_quirk_record_sequence_provenance_and_unpolished_tsv_roundtrip(tmp_path: Path):
    """Verify QuirkRecord sequence_provenance and unpolished fields in TSV SerDe."""
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
            representative_sequence_selected=True,
            sequence_provenance="ONT-derived",
            unpolished=True,
            notes="Bridged variable satellite gap",
        ),
        QuirkRecord(
            feature_id="low_support_stub:V:3p:1",
            chrom="V",
            side="3p",
            quirk_type=IntervalQuirk.LOW_SUPPORT_STUB,
            unitig_ids=("utg179",),
            read_count=2,
            representative_sequence_selected=False,
            sequence_provenance=None,
            unpolished=False,
            notes="Unbridged dropout stub",
        ),
    ]

    buf = StringIO()
    emit_quirks_tsv(records, buf)
    buf.seek(0)
    parsed = parse_quirks_tsv(buf)

    assert len(parsed) == 2
    assert parsed[0].feature_id == "sat_gap_variable:IV:3p:1"
    assert parsed[0].representative_sequence_selected is True
    assert parsed[0].sequence_provenance == "ONT-derived"
    assert parsed[0].unpolished is True

    assert parsed[1].feature_id == "low_support_stub:V:3p:1"
    assert parsed[1].representative_sequence_selected is False
    assert parsed[1].sequence_provenance is None
    assert parsed[1].unpolished is False


def test_resolver_satellite_gap_and_dropout_populate_ont_provenance_and_unpolished():
    """Verify resolver populates sequence_provenance='ONT-derived' and unpolished=True when rep selected."""
    engine = OntBridgingEngine()

    p_left = SegmentPort("utg_l", Strand.FORWARD)
    p_right = SegmentPort("utg_r", Strand.FORWARD)

    # 1. resolve_satellite_gap with bridging reads
    paf_lines = [
        "read1\t30000\t1000\t3000\t+\tutg_l\t50000\t48000\t50000\t2000\t2000\t60",
        "read1\t30000\t4000\t6000\t+\tutg_r\t50000\t0\t2000\t2000\t2000\t60",
        "read2\t30000\t1000\t3000\t+\tutg_l\t50000\t48000\t50000\t2000\t2000\t60",
        "read2\t30000\t4000\t6000\t+\tutg_r\t50000\t0\t2000\t2000\t2000\t60",
    ]
    reads_seq = {
        "read1": "A" * 3000 + "G" * 1000 + "T" * 26000,
        "read2": "A" * 3000 + "G" * 1000 + "T" * 26000,
    }
    sat_res = engine.resolve_satellite_gap(
        chrom="IV",
        side="3p",
        left_port=p_left,
        right_port=p_right,
        paf_source=paf_lines,
        reads_sequences=reads_seq,
    )
    assert sat_res.quirk_record.representative_sequence_selected is True
    assert sat_res.quirk_record.sequence_provenance == "ONT-derived"
    assert sat_res.quirk_record.unpolished is True

    # 2. resolve_somatic_dropout with bridging reads
    p_som = SegmentPort("utg_som", Strand.REVERSE)
    drop_quirk, drop_end, rep_seq = engine.resolve_somatic_dropout(
        chrom="IV",
        side="3p",
        stem_port=p_left,
        germline_port=p_right,
        somatic_port=p_som,
        paf_path=paf_lines,
        reads_fasta=reads_seq,
    )
    assert drop_quirk.representative_sequence_selected is True
    assert drop_quirk.sequence_provenance == "ONT-derived"
    assert drop_quirk.unpolished is True
    assert rep_seq is not None


def test_fasta_agp_emitter_unknown_gap_rdna_tangle():
    """Verify ScaffoldStep with is_unknown_gap=True emits 100 'N's in FASTA and U 100 repeat yes rna_tangle in AGP."""
    emitter = FastaAgpEmitter()
    step = ScaffoldStep(
        segment_id="",
        strand=Strand.FORWARD,
        provenance=SequenceProvenance.RDNA_CHAINED,
        is_unknown_gap=True,
    )
    plan = ScaffoldPlan(scaffold_id="chrII", steps=[step])
    seq, agp_rows = emitter.emit_scaffold(plan)

    assert seq == "N" * 100
    assert len(agp_rows) == 1

    row = agp_rows[0]
    assert isinstance(row, AgpGapRow)
    assert row.object_id == "chrII"
    assert row.object_beg == 1
    assert row.object_end == 100
    assert row.part_number == 1
    assert row.component_type == "U"
    assert row.gap_length == 100
    assert row.gap_type == "repeat"
    assert row.linkage == "yes"
    assert row.linkage_evidence == "rna_tangle"

    line = row.to_line()
    assert line == "chrII\t1\t100\t1\tU\t100\trepeat\tyes\trna_tangle"

    # Verify convenience function emit_scaffold
    conv_seq, conv_rows = emit_scaffold(plan)
    assert conv_seq == "N" * 100
    assert len(conv_rows) == 1

    buf_fa = StringIO()
    buf_agp = StringIO()
    emitter.emit_scaffold(plan, fasta_output=buf_fa, agp_output=buf_agp)
    assert buf_fa.getvalue() == ">chrII\n" + ("N" * 100) + "\n"
    assert "chrII\t1\t100\t1\tU\t100\trepeat\tyes\trna_tangle\n" in buf_agp.getvalue()

    # Multi-step plan with flanking unitigs
    hifi_step = ScaffoldStep(
        segment_id="stem_01",
        strand="+",
        sequence="A" * 500,
        provenance=SequenceProvenance.HIFI_CORE,
    )
    cap_step = ScaffoldStep(
        segment_id="cap_01",
        strand="+",
        sequence="T" * 300,
        provenance=SequenceProvenance.RDNA_CHAINED,
    )
    full_plan = ScaffoldPlan(scaffold_id="chrII", steps=[hifi_step, step, cap_step])
    full_seq, full_rows = emitter.emit_scaffold(full_plan)
    assert full_seq == ("A" * 500) + ("N" * 100) + ("T" * 300)
    assert len(full_rows) == 3
    assert full_rows[1].to_line() == "chrII\t501\t600\t2\tU\t100\trepeat\tyes\trna_tangle"


def test_path_component_default_hifi_spanning_and_gap_emission():
    """Verify PathComponent defaults to hifi-spanning linkage evidence and emits U-gap in AGP."""
    comp = PathComponent(segment_id="utg01", strand="+", sequence="ACGT")
    assert comp.linkage_evidence == "hifi-spanning"

    c1 = PathComponent(segment_id="s1", strand="+", sequence="A" * 200, provenance=SequenceProvenance.HIFI_CORE)
    gap_comp = PathComponent(
        is_gap=True,
        gap_length=100,
        gap_type="repeat",
        linkage="yes",
        linkage_evidence="hifi-spanning",
        gap_component_type="N",
    )
    c2 = PathComponent(segment_id="s2", strand="+", sequence="C" * 200, provenance=SequenceProvenance.HIFI_CORE)
    seq, rows = assemble_scaffold_sequence([c1, gap_comp, c2], scaffold_id="chrX")
    assert len(rows) == 3
    assert isinstance(rows[1], AgpGapRow)
    assert rows[1].component_type == "U"
    assert rows[1].linkage_evidence == "hifi-spanning"
    assert rows[1].to_line() == "chrX\t201\t300\t2\tU\t100\trepeat\tyes\thifi-spanning"


def test_emit_grs_bed_coordinates_and_strict_exclusions(tmp_path: Path):
    """Verify emit_grs_bed emits 0-based half-open coordinates and strictly excludes repeat gaps and ONT."""
    recs = [
        # HiFi-core segment (should be excluded)
        AgpRecord(
            object_id="chrI",
            object_beg=1,
            object_end=1000,
            part_number=1,
            component_type="W",
            component_id="utg_core",
            component_beg=1,
            component_end=1000,
            orientation="+",
            sequence_provenance=SequenceProvenance.HIFI_CORE,
        ),
        # HiFi-germline GRS segment (should be included)
        AgpRecord(
            object_id="chrI",
            object_beg=1001,
            object_end=2500,
            part_number=2,
            component_type="W",
            component_id="utg_grs",
            component_beg=1,
            component_end=1500,
            orientation="-",
            sequence_provenance=SequenceProvenance.HIFI_GERMLINE,
        ),
        # Gap type repeat / variable satellite gap (strictly excluded)
        AgpGapRow(
            object_id="chrI",
            object_beg=2501,
            object_end=2600,
            part_number=3,
            component_type="U",
            gap_length=100,
            gap_type="repeat",
            linkage="yes",
            linkage_evidence="satellite-kmer",
        ),
        # ONT-derived unpolished bridge (strictly excluded)
        AgpRecord(
            object_id="chrI",
            object_beg=2601,
            object_end=3600,
            part_number=4,
            component_type="W",
            component_id="ont_read_01",
            component_beg=1,
            component_end=1000,
            orientation="+",
            sequence_provenance=SequenceProvenance.ONT_DERIVED,
            unpolished=True,
        ),
        # Second HiFi-germline segment on chrII
        AgpRecord(
            object_id="chrII",
            object_beg=5001,
            object_end=7000,
            part_number=1,
            component_type="W",
            component_id="utg_grs2",
            component_beg=1,
            component_end=2000,
            orientation="+",
            sequence_provenance=SequenceProvenance.HIFI_GERMLINE,
        ),
    ]

    out_bed = tmp_path / "GRS.bed"
    ret = emit_grs_bed(recs, out_bed)
    assert ret == out_bed
    assert out_bed.exists()

    lines = [line.strip().split("\t") for line in out_bed.read_text(encoding="utf-8").splitlines() if line.strip()]
    # Only the two HiFi-germline segments should be present
    assert len(lines) == 2

    # Verify 0-based half-open coordinates:
    # chrI: 1001..2500 in AGP -> 1000..2500 in BED
    assert lines[0] == ["chrI", "1000", "2500", "utg_grs", "0", "-"]
    # chrII: 5001..7000 in AGP -> 5000..7000 in BED
    assert lines[1] == ["chrII", "5000", "7000", "utg_grs2", "0", "+"]


def test_emit_telomere_bed_5p_3p_caps(tmp_path: Path):
    """Verify emit_telomere_bed accurately detects and annotates 5' and 3' terminal telomere caps in BED format."""
    motif_fwd = DEFAULT_POLICY.telomere.motif_forward.value  # TTAGGC
    motif_rev = DEFAULT_POLICY.telomere.motif_reverse.value  # GCCTAA

    # chr1: 5' cap (GCCTAA*35 at offset 15) and 3' cap (TTAGGC*35 at suffix offset 0)
    seq_chr1 = ("A" * 15) + (motif_rev * 35) + ("C" * 2000) + (motif_fwd * 35)
    # chr2: 5' cap with fallback forward motif (TTAGGC*35 at offset 0)
    seq_chr2 = (motif_fwd * 35) + ("G" * 3000)
    # chr3: naked sequence without terminal telomeres
    seq_chr3 = "A" * 4000

    scaffolds = {
        "chr1": seq_chr1,
        "chr2": seq_chr2,
        "chr3": seq_chr3,
    }

    out_bed = tmp_path / "telomeres.bed"
    ret = emit_telomere_bed(scaffolds, out_bed)
    assert ret == out_bed
    assert out_bed.exists()

    lines = [line.strip().split("\t") for line in out_bed.read_text(encoding="utf-8").splitlines() if line.strip()]
    assert len(lines) == 3

    # chr1 5p cap: offset 15, length 35*6 = 210 -> chromStart=15, chromEnd=225, strand="-"
    assert lines[0] == ["chr1", "15", "225", "telomere_5p", "35", "-"]
    # chr1 3p cap: ends at len(seq_chr1), length 210 -> chromStart=len-210, chromEnd=len, strand="+"
    chr1_len = len(seq_chr1)
    assert lines[1] == ["chr1", str(chr1_len - 210), str(chr1_len), "telomere_3p", "35", "+"]
    # chr2 5p cap: offset 0, fallback motif TTAGGC -> strand="+"
    assert lines[2] == ["chr2", "0", "210", "telomere_5p", "35", "+"]


def test_find_break_site_scaffolds_and_nominate_super():
    """Verify nomination of SUPER candidate chromosomes sorted descending by length."""
    scaffolds = {
        "chr_big": "A" * 500000,
        "chr_mid": "C" * 300000,
        "chr_small": "G" * 50000,  # Below 100 kb threshold
        "chr_bacterial": "T" * 600000,  # No break site
    }

    agp_records = [
        AgpRecord(
            object_id="chr_big",
            object_beg=1,
            object_end=500000,
            part_number=1,
            component_type="W",
            component_id="utg_big",
            component_beg=1,
            component_end=500000,
            orientation="+",
            sequence_provenance=SequenceProvenance.HIFI_CORE,
        ),
        AgpRecord(
            object_id="chr_mid",
            object_beg=1,
            object_end=300000,
            part_number=1,
            component_type="W",
            component_id="utg_mid",
            component_beg=1,
            component_end=300000,
            orientation="+",
            sequence_provenance=SequenceProvenance.HIFI_CORE,
        ),
        AgpRecord(
            object_id="chr_small",
            object_beg=1,
            object_end=50000,
            part_number=1,
            component_type="W",
            component_id="utg_small",
            component_beg=1,
            component_end=50000,
            orientation="+",
            sequence_provenance=SequenceProvenance.HIFI_CORE,
        ),
        AgpRecord(
            object_id="chr_bacterial",
            object_beg=1,
            object_end=600000,
            part_number=1,
            component_type="W",
            component_id="utg_bact",
            component_beg=1,
            component_end=600000,
            orientation="+",
            sequence_provenance=SequenceProvenance.HIFI_CORE,
        ),
    ]

    junctions = [
        JunctionRecord(
            feature_id="junc_big",
            chrom="chr_big",
            side="3p",
            junction_class=JunctionClass.PDE_TERMINAL,
            stem_port="utg_big+",
            germline_continuation_port="utg_big_grs+",
        ),
        JunctionRecord(
            feature_id="junc_mid",
            chrom="chr_mid",
            side="5p",
            junction_class=JunctionClass.PDE_TERMINAL,
            stem_port="utg_mid-",
        ),
        JunctionRecord(
            feature_id="junc_small",
            chrom="chr_small",
            side="3p",
            junction_class=JunctionClass.PDE_TERMINAL,
            stem_port="utg_small+",
        ),
    ]

    resolutions = [
        ResolutionRecord(
            feature_id="res_big",
            chrom="chr_big",
            side="3p",
            resolution_type=ResolutionType.PDE_TERMINAL_CAP,
        ),
        ResolutionRecord(
            feature_id="res_bact",
            chrom="chr_bacterial",
            side="3p",
            resolution_type=ResolutionType.PDE_TERMINAL_CAP,
        ),
    ]

    new_scaffs, new_agp, rename_map = nominate_super_scaffolds(
        scaffolds=scaffolds,
        agp_records=agp_records,
        junctions=junctions,
        min_super_length=100000,
        resolutions=resolutions,
    )

    # chr_big (500 kb) qualifies and is largest -> SUPER_1
    # chr_mid (300 kb) qualifies and is second largest -> SUPER_2
    # chr_small (50 kb) < 100 kb -> retains chr_small
    # chr_bacterial (600 kb) has NO break site -> retains chr_bacterial
    assert rename_map == {
        "chr_big": "SUPER_1",
        "chr_mid": "SUPER_2",
    }

    assert "SUPER_1" in new_scaffs and len(new_scaffs["SUPER_1"]) == 500000
    assert "SUPER_2" in new_scaffs and len(new_scaffs["SUPER_2"]) == 300000
    assert "chr_small" in new_scaffs
    assert "chr_bacterial" in new_scaffs

    # AGP records updated
    obj_ids = [r.object_id for r in new_agp]
    assert obj_ids == ["SUPER_1", "SUPER_2", "chr_small", "chr_bacterial"]

    # Resolution records updated
    assert resolutions[0].chrom == "SUPER_1"
    assert resolutions[1].chrom == "chr_bacterial"


def test_compute_break_site_summary_and_emit_summary_tsv(tmp_path: Path):
    """Verify compute_break_site_summary and emit_summary_tsv compute accurate core vs GRS stoichiometry."""
    graph = BidirectedGraph()
    # chr1 components: 100 kb core (60x), 30 kb GRS (6x)
    graph.add_segment("chr1_stem", length=100000, depth=60.0)
    graph.add_segment("chr1_grs", length=30000, depth=6.0)
    # chr2 components: 200 kb core (50x), 20 kb GRS (10x)
    graph.add_segment("chr2_stem", length=200000, depth=50.0)
    graph.add_segment("chr2_grs", length=20000, depth=10.0)
    # Bacterial contaminant (no break site): 500 kb core (15x)
    graph.add_segment("bact_utg", length=500000, depth=15.0)

    scaffolds = {
        "SUPER_1": "A" * 220000,
        "SUPER_2": "C" * 130000,
        "bacterial_contig": "G" * 500000,
    }

    agp_records = [
        # SUPER_1: chr2
        AgpRecord(
            object_id="SUPER_1",
            object_beg=1,
            object_end=200000,
            part_number=1,
            component_type="W",
            component_id="chr2_stem",
            component_beg=1,
            component_end=200000,
            orientation="+",
            sequence_provenance=SequenceProvenance.HIFI_CORE,
        ),
        AgpRecord(
            object_id="SUPER_1",
            object_beg=200001,
            object_end=220000,
            part_number=2,
            component_type="W",
            component_id="chr2_grs",
            component_beg=1,
            component_end=20000,
            orientation="+",
            sequence_provenance=SequenceProvenance.HIFI_GERMLINE,
        ),
        # SUPER_2: chr1
        AgpRecord(
            object_id="SUPER_2",
            object_beg=1,
            object_end=100000,
            part_number=1,
            component_type="W",
            component_id="chr1_stem",
            component_beg=1,
            component_end=100000,
            orientation="+",
            sequence_provenance=SequenceProvenance.HIFI_CORE,
        ),
        AgpRecord(
            object_id="SUPER_2",
            object_beg=100001,
            object_end=130000,
            part_number=2,
            component_type="W",
            component_id="chr1_grs",
            component_beg=1,
            component_end=30000,
            orientation="+",
            sequence_provenance=SequenceProvenance.HIFI_GERMLINE,
        ),
        # Bacterial contig
        AgpRecord(
            object_id="bacterial_contig",
            object_beg=1,
            object_end=500000,
            part_number=1,
            component_type="W",
            component_id="bact_utg",
            component_beg=1,
            component_end=500000,
            orientation="+",
            sequence_provenance=SequenceProvenance.HIFI_CORE,
        ),
    ]

    junctions = [
        JunctionRecord(
            feature_id="j1",
            chrom="SUPER_2",
            side="3p",
            junction_class=JunctionClass.PDE_TERMINAL,
            stem_port="chr1_stem+",
            germline_continuation_port="chr1_grs+",
        ),
        JunctionRecord(
            feature_id="j2",
            chrom="SUPER_1",
            side="3p",
            junction_class=JunctionClass.PDE_TERMINAL,
            stem_port="chr2_stem+",
            germline_continuation_port="chr2_grs+",
        ),
    ]

    summary = compute_break_site_summary(graph, scaffolds, agp_records, junctions)

    # Bacterial contig must be strictly excluded!
    assert len(summary.qualifying_sequences) == 2
    s1, s2 = summary.qualifying_sequences
    assert s1.sequence_id == "SUPER_1"
    assert s1.core_length_bp == 200000
    assert s1.core_depth == 50.0
    assert s1.grs_length_bp == 20000
    assert s1.grs_depth == 10.0
    assert s1.stoichiometry_ratio == pytest.approx(10.0 / 50.0)
    assert s1.grs_percentage == pytest.approx(20000 / 220000 * 100.0)
    assert s1.break_site_count == 1

    assert s2.sequence_id == "SUPER_2"
    assert s2.core_length_bp == 100000
    assert s2.core_depth == 60.0
    assert s2.grs_length_bp == 30000
    assert s2.grs_depth == 6.0
    assert s2.stoichiometry_ratio == pytest.approx(6.0 / 60.0)
    assert s2.grs_percentage == pytest.approx(30000 / 130000 * 100.0)
    assert s2.break_site_count == 1

    # Genome-wide totals across qualifying sequences:
    # Total Core = 300,000 bp; mean depth = (200000*50 + 100000*60) / 300000 = 16000000 / 300000 = 53.333x
    assert summary.total_core_length_bp == 300000
    assert summary.mean_core_depth == pytest.approx(53.3333, rel=1e-3)
    # Total GRS = 50,000 bp; mean depth = (20000*10 + 30000*6) / 50000 = 380000 / 50000 = 7.60x
    assert summary.total_grs_length_bp == 50000
    assert summary.mean_grs_depth == pytest.approx(7.60, rel=1e-3)
    # Stoichiometry ratio = 7.60 / 53.3333 = 0.1425
    assert summary.overall_stoichiometry_ratio == pytest.approx(7.60 / (160 / 3), rel=1e-3)
    # GRS percentage = 50000 / 350000 * 100 = 14.2857%
    assert summary.overall_grs_percentage == pytest.approx(50000 / 350000 * 100.0, rel=1e-3)
    assert summary.total_break_sites == 2

    # Emission and parsing round-trip
    out_tsv = tmp_path / "summary.tsv"
    emit_summary_tsv(summary, out_tsv)
    assert out_tsv.exists()

    seq_rows, total_row = parse_summary_tsv(out_tsv)
    assert len(seq_rows) == 2
    assert seq_rows[0]["sequence_id"] == "SUPER_1"
    assert seq_rows[1]["sequence_id"] == "SUPER_2"
    assert total_row["sequence_id"] == "TOTAL"
    assert total_row["core_length_bp"] == 300000
    assert total_row["grs_length_bp"] == 50000

    # Human-readable table rendering
    table_str = format_summary_table(summary)
    assert "gfapdefinisher Finish Summary" in table_str
    assert "SUPER_1" in table_str
    assert "SUPER_2" in table_str
    assert "TOTAL" in table_str
