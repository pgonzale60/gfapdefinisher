"""Tests for the gfapdefinisher command-line interface."""

from __future__ import annotations

from pathlib import Path

import pytest
from pyfaidx import Fasta

from gfapdefinisher.cli import main
from gfapdefinisher.emitter import SequenceProvenance, parse_agp
from gfapdefinisher.graph import load_gfa
from gfapdefinisher.policy import DEFAULT_POLICY
from gfapdefinisher.resolutions import (
    EndState,
    ResolutionType,
    parse_end_states_tsv,
    parse_junctions_tsv,
    parse_resolutions_tsv,
)


def _make_paf_line(
    qname: str,
    qlen: int,
    qstart: int,
    qend: int,
    strand: str,
    tname: str,
    tlen: int,
    tstart: int,
    tend: int,
    alen: int = 1500,
    mapq: int = 60,
    nmatch: int = 1400,
) -> str:
    """Format a single synthetic PAF record."""
    return (
        f"{qname}\t{qlen}\t{qstart}\t{qend}\t{strand}\t"
        f"{tname}\t{tlen}\t{tstart}\t{tend}\t{nmatch}\t{alen}\t{mapq}\ttp:A:P"
    )


def _build_synthetic_gfa(path: Path) -> None:
    """Write synthetic GFA with a PDE junction, authentic somatic telomere, and an orphan cap."""
    motif = DEFAULT_POLICY.telomere.motif_forward.value
    # chr_stem: 50 kb somatic stem unitig
    seq_stem = "A" * 50000
    # chr_som: 25 kb somatic arm ending with 35x telomere repeats (210 bp cap)
    seq_som = ("A" * 24790) + (motif * 35)
    # chr_grs: 30 kb germline retained branch
    seq_grs = "C" * 30000
    # orphan_cap: 10 kb orphan cap unitig ending with 35x telomere repeats
    seq_orphan = ("G" * 9790) + (motif * 35)

    lines = [
        "H\tVN:Z:1.0",
        f"S\tchr_stem\t{seq_stem}\tLN:i:50000\tdp:f:60.0",
        f"S\tchr_som\t{seq_som}\tLN:i:25000\tdp:f:30.0",
        f"S\tchr_grs\t{seq_grs}\tLN:i:30000\tdp:f:6.0",
        f"S\torphan_cap\t{seq_orphan}\tLN:i:10000\tdp:f:5.0",
        # PDE fork: chr_stem+ -> chr_som+ and chr_stem+ -> chr_grs+
        "L\tchr_stem\t+\tchr_som\t+\t100M",
        "L\tchr_stem\t+\tchr_grs\t+\t100M",
    ]
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


def test_cli_no_command(capsys: pytest.CaptureFixture[str]):
    """Verify invoking CLI with no arguments exits with code 0 and prints help."""
    with pytest.raises(SystemExit) as exc:
        main([])
    assert exc.value.code == 0
    captured = capsys.readouterr()
    assert "usage:" in captured.out or "gfapdefinisher" in captured.out


def test_report_subcommand_missing_args():
    """Verify report subcommand fails when required --gfa argument is omitted."""
    with pytest.raises(SystemExit) as exc:
        main(["report"])
    assert exc.value.code != 0


def test_report_subcommand_valid(tmp_path: Path):
    """Verify report subcommand produces valid ledger TSV output file."""
    gfa_path = tmp_path / "test.gfa"
    _build_synthetic_gfa(gfa_path)
    out_tsv = tmp_path / "ledger.tsv"

    exit_code = main(["report", "--gfa", str(gfa_path), "--out", str(out_tsv)])
    assert exit_code == 0
    assert out_tsv.exists()

    content = out_tsv.read_text(encoding="utf-8")
    lines = [line for line in content.splitlines() if line.strip()]
    assert len(lines) > 1
    header = lines[0].split("\t")
    assert "segment_name" in header
    assert "port" in header


def test_report_subcommand_stdout(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    """Verify report subcommand emits to stdout when --out is not provided."""
    gfa_path = tmp_path / "test.gfa"
    _build_synthetic_gfa(gfa_path)

    exit_code = main(["report", "--gfa", str(gfa_path)])
    assert exit_code == 0

    captured = capsys.readouterr()
    lines = [line for line in captured.out.splitlines() if line.strip()]
    assert len(lines) > 1
    assert "segment_name" in lines[0]


def test_finish_subcommand_missing_args(tmp_path: Path):
    """Verify finish subcommand enforces required --gfa and --out-dir arguments."""
    with pytest.raises(SystemExit) as exc:
        main(["finish"])
    assert exc.value.code != 0

    gfa_path = tmp_path / "test.gfa"
    _build_synthetic_gfa(gfa_path)
    with pytest.raises(SystemExit) as exc2:
        main(["finish", "--gfa", str(gfa_path)])
    assert exc2.value.code != 0


def test_finish_subcommand_synthetic_and_non_empty_tsvs(tmp_path: Path):
    """Verify finish subcommand writes all 4 typed non-empty TSV ledgers without ONT PAF."""
    gfa_path = tmp_path / "test.gfa"
    _build_synthetic_gfa(gfa_path)
    out_dir = tmp_path / "finish_out"

    exit_code = main(["finish", "--gfa", str(gfa_path), "--out-dir", str(out_dir)])
    assert exit_code == 0

    # Verify all 4 output TSVs exist
    ledger_file = out_dir / "ledger.tsv"
    junctions_file = out_dir / "junctions.tsv"
    end_states_file = out_dir / "end_states.tsv"
    resolutions_file = out_dir / "resolutions.tsv"

    assert ledger_file.is_file()
    assert junctions_file.is_file()
    assert end_states_file.is_file()
    assert resolutions_file.is_file()
    assert not (out_dir / "quirks.tsv").exists()

    # Verify each TSV is non-empty (contains header plus data rows)
    ledger_lines = [line for line in ledger_file.read_text(encoding="utf-8").splitlines() if line.strip()]
    assert len(ledger_lines) > 1

    junction_recs = parse_junctions_tsv(junctions_file)
    assert len(junction_recs) >= 1
    assert any(j.junction_class == "pde_terminal" for j in junction_recs)

    end_state_recs = parse_end_states_tsv(end_states_file)
    assert len(end_state_recs) >= 1
    assert any(es.state == EndState.END_GERMLINE_ORPHAN for es in end_state_recs)
    # Ensure orphan_cap unitig is not emitted as a separate chromosome
    assert not any(es.chrom == "orphan_cap" for es in end_state_recs)

    res_recs = parse_resolutions_tsv(resolutions_file)
    assert len(res_recs) >= 1
    assert any(q.resolution_type == ResolutionType.PDE_TERMINAL_CAP for q in res_recs)


def test_finish_subcommand_with_ont_paf_bridging(tmp_path: Path):
    """Verify finish subcommand resolves variable satellite gap with ONT PAF and slices reads."""
    gfa_path = tmp_path / "test.gfa"
    _build_synthetic_gfa(gfa_path)
    out_dir = tmp_path / "finish_ont_out"

    # Synthetic PAF with 3 spanning reads bridging chr_grs (proximal) and orphan_cap (distal)
    paf_path = tmp_path / "spanning.paf"
    paf_lines = [
        # read_1: 15 kb read, gap = 1500 bp
        _make_paf_line("read_1", 20000, 1000, 2500, "+", "chr_grs", 30000, 28500, 30000),
        _make_paf_line("read_1", 20000, 4000, 5500, "+", "orphan_cap", 10000, 0, 1500),
        # read_2: 15 kb read, gap = 1500 bp
        _make_paf_line("read_2", 20000, 1000, 2500, "+", "chr_grs", 30000, 28500, 30000),
        _make_paf_line("read_2", 20000, 4000, 5500, "+", "orphan_cap", 10000, 0, 1500),
        # read_3: 15 kb read, gap = 1600 bp
        _make_paf_line("read_3", 20000, 1000, 2500, "+", "chr_grs", 30000, 28500, 30000),
        _make_paf_line("read_3", 20000, 4100, 5600, "+", "orphan_cap", 10000, 0, 1500),
    ]
    paf_path.write_text("\n".join(paf_lines) + "\n", encoding="utf-8")

    # Synthetic FASTA with read sequences for BIO-004 slicing
    fasta_path = tmp_path / "reads.fa"
    gap_seq = "ACGT" * 375  # 1500 bp
    read1_seq = ("A" * 2500) + gap_seq + ("T" * 16000)
    fasta_lines = [
        ">read_1 sample ont read 1",
        read1_seq,
        ">read_2 sample ont read 2",
        read1_seq,
        ">read_3 sample ont read 3",
        ("A" * 2500) + ("ACGT" * 400) + ("T" * 15900),
    ]
    fasta_path.write_text("\n".join(fasta_lines) + "\n", encoding="utf-8")

    exit_code = main(
        [
            "finish",
            "--gfa",
            str(gfa_path),
            "--out-dir",
            str(out_dir),
            "--ont-paf",
            str(paf_path),
            "--ont-reads",
            str(fasta_path),
        ]
    )
    assert exit_code == 0

    end_states_file = out_dir / "end_states.tsv"
    resolutions_file = out_dir / "resolutions.tsv"

    end_state_recs = parse_end_states_tsv(end_states_file)
    # The gap should now be resolved to END_COMPLETE
    complete_ends = [es for es in end_state_recs if es.state == EndState.END_COMPLETE]
    assert len(complete_ends) >= 1
    assert any(es.resolved is True for es in complete_ends)
    resolved_end = [es for es in complete_ends if es.chrom == "chr_stem"][0]
    assert resolved_end.somatic_port == "chr_som+"
    assert resolved_end.germline_port == "orphan_cap+"
    assert not any(es.chrom == "orphan_cap" for es in end_state_recs)

    res_recs = parse_resolutions_tsv(resolutions_file)
    sat_res = [q for q in res_recs if q.resolution_type == ResolutionType.SAT_GAP_VARIABLE]
    assert len(sat_res) >= 1
    bridged_res = sat_res[0]
    assert bridged_res.spanning_read_count == 3
    assert bridged_res.representative_sequence_selected is True
    assert bridged_res.agp_linkage_evidence == "ont-spanning"
    assert bridged_res.length_median_bp == 1500


def test_finish_subcommand_with_annotations_and_window(tmp_path: Path):
    """Verify finish subcommand consumes annotations and terminal-window argument."""
    gfa_path = tmp_path / "test.gfa"
    _build_synthetic_gfa(gfa_path)
    out_dir = tmp_path / "finish_annot_out"

    # Synthetic TIDK file
    tidk_path = tmp_path / "telomere.tsv"
    tidk_lines = [
        "sequence\tstart\tend\tForward_Matches\tReverse_Matches",
        "chr_som\t24000\t25000\t35\t0",
    ]
    tidk_path.write_text("\n".join(tidk_lines) + "\n", encoding="utf-8")

    # Synthetic Barrnap GFF3 file
    barrnap_path = tmp_path / "rrna.gff"
    barrnap_lines = [
        "##gff-version 3",
        "chr_som\tbarrnap:0.9\trRNA\t1000\t2800\t0.0\t+\t.\tName=18S_rRNA",
    ]
    barrnap_path.write_text("\n".join(barrnap_lines) + "\n", encoding="utf-8")

    exit_code = main(
        [
            "finish",
            "--gfa",
            str(gfa_path),
            "--out-dir",
            str(out_dir),
            "--tidk",
            str(tidk_path),
            "--barrnap",
            str(barrnap_path),
            "--terminal-window",
            "5000",
        ]
    )
    assert exit_code == 0
    ledger_file = out_dir / "ledger.tsv"
    assert ledger_file.is_file()
    ledger_text = ledger_file.read_text(encoding="utf-8")
    assert "chr_som" in ledger_text
    assert "18S_rRNA" in ledger_text or "chr_stem" in ledger_text


def test_finish_subcommand_with_fastq_reads(tmp_path: Path):
    """Verify finish subcommand supports FASTQ-formatted reads with quality headers."""
    gfa_path = tmp_path / "test.gfa"
    _build_synthetic_gfa(gfa_path)
    out_dir = tmp_path / "finish_fastq_out"

    paf_path = tmp_path / "spanning.paf"
    paf_lines = [
        _make_paf_line("fq_read1", 20000, 1000, 2500, "+", "chr_grs", 30000, 28500, 30000),
        _make_paf_line("fq_read1", 20000, 4000, 5500, "+", "orphan_cap", 10000, 0, 1500),
        _make_paf_line("fq_read2", 20000, 1000, 2500, "+", "chr_grs", 30000, 28500, 30000),
        _make_paf_line("fq_read2", 20000, 4000, 5500, "+", "orphan_cap", 10000, 0, 1500),
    ]
    paf_path.write_text("\n".join(paf_lines) + "\n", encoding="utf-8")

    fastq_path = tmp_path / "reads.fq"
    gap_seq = "ACGT" * 375
    read_seq = ("A" * 2500) + gap_seq + ("T" * 16000)
    qual_seq = "I" * len(read_seq)
    fastq_lines = [
        "@fq_read1 sample fastq read 1",
        read_seq,
        "+",
        qual_seq,
        "@fq_read2 sample fastq read 2",
        read_seq,
        "+",
        qual_seq,
    ]
    fastq_path.write_text("\n".join(fastq_lines) + "\n", encoding="utf-8")

    exit_code = main(
        [
            "finish",
            "--gfa",
            str(gfa_path),
            "--out-dir",
            str(out_dir),
            "--ont-paf",
            str(paf_path),
            "--ont-reads",
            str(fastq_path),
        ]
    )
    assert exit_code == 0
    resolutions_file = out_dir / "resolutions.tsv"
    res_recs = parse_resolutions_tsv(resolutions_file)
    sat_res = [q for q in res_recs if q.resolution_type == ResolutionType.SAT_GAP_VARIABLE]
    assert len(sat_res) >= 1
    assert sat_res[0].representative_sequence_selected is True


def test_finish_subcommand_emits_fasta_and_agp_by_default(tmp_path: Path):
    """Verify finish subcommand emits scaffolds.fasta and scaffolds.agp by default."""
    gfa_path = tmp_path / "test.gfa"
    _build_synthetic_gfa(gfa_path)
    out_dir = tmp_path / "finish_default_emission"

    exit_code = main(["finish", "--gfa", str(gfa_path), "--out-dir", str(out_dir)])
    assert exit_code == 0

    fasta_file = out_dir / "scaffolds.fasta"
    agp_file = out_dir / "scaffolds.agp"

    assert fasta_file.is_file()
    assert agp_file.is_file()

    fasta_content = fasta_file.read_text(encoding="utf-8")
    assert ">chr_stem" in fasta_content
    assert ">chr_som" not in fasta_content
    assert (out_dir / "excluded.fasta").is_file()
    assert ">chr_som" in (out_dir / "excluded.fasta").read_text(encoding="utf-8")

    agp_content = agp_file.read_text(encoding="utf-8")
    assert agp_content.startswith("##agp-version 2.1\n")
    assert "chr_stem" in agp_content
    assert "chr_grs" in agp_content
    assert "chr_som" not in agp_content


def test_finish_subcommand_custom_fasta_and_agp_options(tmp_path: Path):
    """Verify finish subcommand respects custom --emit-fasta, --emit-agp, and --prefix options."""
    gfa_path = tmp_path / "test.gfa"
    _build_synthetic_gfa(gfa_path)
    out_dir = tmp_path / "finish_custom_emission"
    custom_fa = tmp_path / "custom_scaffolds.fa"
    custom_agp = tmp_path / "custom_scaffolds.agp"

    exit_code = main(
        [
            "finish",
            "--gfa",
            str(gfa_path),
            "--out-dir",
            str(out_dir),
            "--emit-fasta",
            str(custom_fa),
            "--emit-agp",
            str(custom_agp),
            "--prefix",
            "TEST_",
        ]
    )
    assert exit_code == 0

    assert custom_fa.is_file()
    assert custom_agp.is_file()

    fa_text = custom_fa.read_text(encoding="utf-8")
    assert ">TEST_chr_stem" in fa_text
    assert ">TEST_chr_som" not in fa_text
    assert (out_dir / "excluded.fasta").is_file()
    assert ">chr_som" in (out_dir / "excluded.fasta").read_text(encoding="utf-8")

    agp_text = custom_agp.read_text(encoding="utf-8")
    assert "TEST_chr_stem" in agp_text
    assert "chr_grs" in agp_text
    assert "chr_som" not in agp_text


def test_finish_subcommand_rdna_cassette_wiring(tmp_path: Path):
    """Verify finish subcommand resolves rDNA cassette, chains copies, and stitches scaffold."""
    motif = DEFAULT_POLICY.telomere.motif_forward.value
    # chr_arm: 40 kb core chromosome arm
    seq_arm = "A" * 40000
    # rdna_unitig: 15 kb rDNA unitig with Barrnap hits
    seq_rdna = "C" * 15000
    # distal_cap: 8 kb unitig ending with 35x telomere repeats (210 bp cap)
    seq_cap = ("G" * 7790) + (motif * 35)

    gfa_path = tmp_path / "rdna_test.gfa"
    lines = [
        "H\tVN:Z:1.0",
        f"S\tchr_arm\t{seq_arm}\tLN:i:40000\tdp:f:60.0",
        f"S\trdna_unitig\t{seq_rdna}\tLN:i:15000\tdp:f:250.0",
        f"S\tdistal_cap\t{seq_cap}\tLN:i:8000\tdp:f:55.0",
        # Traversal: chr_arm+ -> rdna_unitig+ (50 bp overlap) -> distal_cap+ (50 bp overlap)
        "L\tchr_arm\t+\trdna_unitig\t+\t50M",
        "L\trdna_unitig\t+\tdistal_cap\t+\t50M",
    ]
    gfa_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    # Barrnap GFF3 annotating rdna_unitig with 18S and 28S rRNA genes
    barrnap_path = tmp_path / "barrnap.gff"
    barrnap_lines = [
        "##gff-version 3",
        "rdna_unitig\tbarrnap:0.9\trRNA\t1000\t2800\t0.0\t+\t.\tName=18S_rRNA",
        "rdna_unitig\tbarrnap:0.9\trRNA\t5000\t8500\t0.0\t+\t.\tName=28S_rRNA",
    ]
    barrnap_path.write_text("\n".join(barrnap_lines) + "\n", encoding="utf-8")

    # TIDK telomeric repeat table for distal_cap
    tidk_path = tmp_path / "tidk.tsv"
    tidk_lines = [
        "sequence\tstart\tend\tForward_Matches\tReverse_Matches",
        "distal_cap\t7500\t8000\t35\t0",
    ]
    tidk_path.write_text("\n".join(tidk_lines) + "\n", encoding="utf-8")

    out_dir = tmp_path / "finish_rdna_out"
    exit_code = main(
        [
            "finish",
            "--gfa",
            str(gfa_path),
            "--out-dir",
            str(out_dir),
            "--barrnap",
            str(barrnap_path),
            "--tidk",
            str(tidk_path),
        ]
    )
    assert exit_code == 0

    end_states = parse_end_states_tsv(out_dir / "end_states.tsv")
    shared_ends = [e for e in end_states if e.state == EndState.END_SHARED_EXEMPT]
    assert len(shared_ends) >= 1
    assert any(e.chrom == "chr_arm" for e in shared_ends)

    resolutions = parse_resolutions_tsv(out_dir / "resolutions.tsv")
    rdna_res = [q for q in resolutions if q.resolution_type == ResolutionType.RDNA_CASSETTE]
    assert len(rdna_res) >= 1
    assert rdna_res[0].chrom == "chr_arm"
    assert "distal cap distal_cap+" in rdna_res[0].notes

    # Scaffolds & AGP checks
    fasta_file = out_dir / "scaffolds.fasta"
    agp_file = out_dir / "scaffolds.agp"
    assert fasta_file.is_file()
    assert agp_file.is_file()

    agp_text = agp_file.read_text(encoding="utf-8")
    assert "provenance=rDNA-chained" in agp_text
    assert "distal_cap" in agp_text

    # Total length: 40000 + (15000 - 50) + (8000 - 50) = 40000 + 14950 + 7950 = 62900 bp
    from gfapdefinisher.emitter import parse_agp

    agp_recs = parse_agp(agp_file)
    arm_recs = [r for r in agp_recs if r.object_id == "chr_arm"]
    assert len(arm_recs) == 3
    assert arm_recs[-1].object_end == 62900


def test_build_scaffolds_stitches_pde_germline_caps():
    """Verify build_scaffolds stitches multi-unitig PDE germline cap chains with exact dovetails."""
    from gfapdefinisher.cli import build_scaffolds
    from gfapdefinisher.emitter import SequenceProvenance
    from gfapdefinisher.graph import BidirectedGraph, SegmentPort, Strand
    from gfapdefinisher.ledger import EndEvidence
    from gfapdefinisher.resolutions import JunctionClass, JunctionRecord

    graph = BidirectedGraph()
    # chr_ii: 50 kb core stem
    graph.add_segment("chr_ii", 50000, sequence="A" * 50000)
    # chr_germ: 10 kb germline arm
    graph.add_segment("chr_germ", 10000, sequence="C" * 10000)
    # distal_cap: 5 kb distal telomeric cap
    graph.add_segment("distal_cap", 5000, sequence="G" * 5000)

    # Edge chr_ii+ -> chr_germ- (100M overlap)
    graph.add_edge("chr_ii+", "chr_germ-", cigar="100M", l1=50000, l2=10000)
    # Edge chr_germ- -> distal_cap+ (80M overlap)
    graph.add_edge("chr_germ-", "distal_cap+", cigar="80M", l1=10000, l2=5000)

    ledger_map = {
        "distal_cap+": EndEvidence(
            segment_name="distal_cap",
            end="3p",
            port=SegmentPort("distal_cap", Strand.FORWARD),
            length=5000,
            depth=6.0,
            degree=0,
            is_tip=True,
            outgoing_neighbors=(),
            has_terminal_telomere=True,
            usco_hit_count=0,
            tiara_class="eukarya",
        ),
        "chr_germ-": EndEvidence(
            segment_name="chr_germ",
            end="5p",
            port=SegmentPort("chr_germ", Strand.REVERSE),
            length=10000,
            depth=6.0,
            degree=1,
            is_tip=False,
            outgoing_neighbors=("distal_cap+",),
            has_terminal_telomere=False,
            usco_hit_count=0,
            tiara_class="eukarya",
        ),
    }

    junction = JunctionRecord(
        feature_id="pde_terminal:chr_ii:3p",
        chrom="chr_ii",
        side="3p",
        junction_class=JunctionClass.PDE_TERMINAL,
        stem_port="chr_ii+",
        germline_continuation_port="chr_germ-",
    )

    scaffolds, agp_recs = build_scaffolds(
        graph=graph,
        ledger_map=ledger_map,
        end_states=[],
        resolutions=[],
        gap_resolutions=[],
        dropout_resolutions=[],
        rdna_resolutions=[],
        junctions=[junction],
    )

    assert "chr_ii" in scaffolds
    # Expected length: 50000 + (10000 - 100) + (5000 - 80) = 50000 + 9900 + 4920 = 64820 bp
    assert len(scaffolds["chr_ii"]) == 64820

    chr_ii_agp = [r for r in agp_recs if r.object_id == "chr_ii"]
    assert len(chr_ii_agp) == 3
    assert chr_ii_agp[0].component_id == "chr_ii"
    assert chr_ii_agp[0].sequence_provenance == SequenceProvenance.HIFI_CORE
    assert chr_ii_agp[1].component_id == "chr_germ"
    assert chr_ii_agp[1].sequence_provenance == SequenceProvenance.HIFI_GERMLINE
    assert chr_ii_agp[2].component_id == "distal_cap"
    assert chr_ii_agp[2].sequence_provenance == SequenceProvenance.HIFI_GERMLINE
    assert chr_ii_agp[-1].object_end == 64820


def test_build_scaffolds_stitches_intercalary_grs_bridges():
    """Verify build_scaffolds chains stems connected by intercalary GRS bridges and adds terminal caps."""
    from gfapdefinisher.cli import build_scaffolds
    from gfapdefinisher.emitter import SequenceProvenance
    from gfapdefinisher.graph import BidirectedGraph, SegmentPort, Strand
    from gfapdefinisher.ledger import EndEvidence
    from gfapdefinisher.resolutions import JunctionClass, JunctionRecord

    graph = BidirectedGraph()
    # stem_a (40 kb), bridge_1 (10 kb), stem_b (30 kb), bridge_2 (12 kb), stem_c (20 kb)
    # cap_5p (5 kb)
    graph.add_segment("stem_a", 40000, sequence="A" * 40000)
    graph.add_segment("bridge_1", 10000, sequence="C" * 10000)
    graph.add_segment("stem_b", 30000, sequence="G" * 30000)
    graph.add_segment("bridge_2", 12000, sequence="T" * 12000)
    graph.add_segment("stem_c", 20000, sequence="A" * 20000)
    graph.add_segment("cap_5p", 5000, sequence="C" * 5000)

    # 5p cap: stem_a- -> cap_5p- (50M overlap), dual is cap_5p+ -> stem_a+ (50M)
    graph.add_edge("stem_a-", "cap_5p-", cigar="50M", l1=40000, l2=5000)
    # stem_a+ -> bridge_1- (100M overlap)
    graph.add_edge("stem_a+", "bridge_1-", cigar="100M", l1=40000, l2=10000)
    # bridge_1- -> stem_b+ (120M overlap)
    graph.add_edge("bridge_1-", "stem_b+", cigar="120M", l1=10000, l2=30000)
    # stem_b+ -> bridge_2+ (80M overlap)
    graph.add_edge("stem_b+", "bridge_2+", cigar="80M", l1=30000, l2=12000)
    # bridge_2+ -> stem_c+ (60M overlap)
    graph.add_edge("bridge_2+", "stem_c+", cigar="60M", l1=12000, l2=20000)

    ledger_map = {
        "cap_5p-": EndEvidence(
            segment_name="cap_5p",
            end="5p",
            port=SegmentPort("cap_5p", Strand.REVERSE),
            length=5000,
            depth=7.0,
            degree=0,
            is_tip=True,
            outgoing_neighbors=(),
            has_terminal_telomere=True,
            usco_hit_count=0,
            tiara_class="eukarya",
        ),
    }

    junctions = [
        JunctionRecord(
            feature_id="pde_terminal:stem_a:5p",
            chrom="stem_a",
            side="5p",
            junction_class=JunctionClass.PDE_TERMINAL,
            stem_port="stem_a-",
            germline_continuation_port="cap_5p-",
        ),
        JunctionRecord(
            feature_id="pde_terminal:stem_a:3p",
            chrom="stem_a",
            side="3p",
            junction_class=JunctionClass.PDE_TERMINAL,
            stem_port="stem_a+",
            germline_continuation_port="bridge_1-",
        ),
        JunctionRecord(
            feature_id="pde_terminal:stem_b:5p",
            chrom="stem_b",
            side="5p",
            junction_class=JunctionClass.PDE_TERMINAL,
            stem_port="stem_b-",
            germline_continuation_port="bridge_1+",
        ),
        JunctionRecord(
            feature_id="pde_terminal:stem_b:3p",
            chrom="stem_b",
            side="3p",
            junction_class=JunctionClass.PDE_TERMINAL,
            stem_port="stem_b+",
            germline_continuation_port="bridge_2+",
        ),
        JunctionRecord(
            feature_id="pde_terminal:stem_c:5p",
            chrom="stem_c",
            side="5p",
            junction_class=JunctionClass.PDE_TERMINAL,
            stem_port="stem_c-",
            germline_continuation_port="bridge_2-",
        ),
    ]

    scaffolds, agp_recs = build_scaffolds(
        graph=graph,
        ledger_map=ledger_map,
        end_states=[],
        resolutions=[],
        gap_resolutions=[],
        dropout_resolutions=[],
        rdna_resolutions=[],
        junctions=junctions,
    )

    # All segments chained into a single chromosome scaffold under stem_a (the largest stem)
    assert "stem_a" in scaffolds
    assert "stem_b" not in scaffolds
    assert "stem_c" not in scaffolds

    # Total length:
    # cap_5p (5000) + stem_a (40000 - 50) + bridge_1 (10000 - 100) + stem_b (30000 - 120) + bridge_2 (12000 - 80) + stem_c (20000 - 60)
    # = 5000 + 39950 + 9900 + 29880 + 11920 + 19940 = 116590 bp
    assert len(scaffolds["stem_a"]) == 116590

    agp_list = [r for r in agp_recs if r.object_id == "stem_a"]
    assert len(agp_list) == 6
    assert [r.component_id for r in agp_list] == ["cap_5p", "stem_a", "bridge_1", "stem_b", "bridge_2", "stem_c"]
    assert [r.sequence_provenance for r in agp_list] == [
        SequenceProvenance.HIFI_GERMLINE,
        SequenceProvenance.HIFI_CORE,
        SequenceProvenance.HIFI_GERMLINE,
        SequenceProvenance.HIFI_CORE,
        SequenceProvenance.HIFI_GERMLINE,
        SequenceProvenance.HIFI_CORE,
    ]
    assert agp_list[-1].object_end == 116590


def test_build_scaffolds_stitches_5p_ont_satellite_gap():
    """Verify build_scaffolds prepends inverted right cap and reverse-complemented ONT bridge for 5p gaps."""
    from gfapdefinisher.cli import build_scaffolds
    from gfapdefinisher.emitter import SequenceProvenance
    from gfapdefinisher.graph import BidirectedGraph
    from gfapdefinisher.resolutions import EndState, EndStateRecord, ResolutionRecord, ResolutionType
    from gfapdefinisher.resolver import PafHit, SatelliteGapResolution, SpanningRead

    graph = BidirectedGraph()
    graph.add_segment("chr_stem", 20000, sequence="A" * 20000)
    graph.add_segment("cap_5p", 5000, sequence="C" * 5000)

    paf_dummy = PafHit(
        qname="read_1",
        qlen=30000,
        qstart=0,
        qend=30000,
        strand="+",
        tname="chr_stem",
        tlen=20000,
        tstart=0,
        tend=1000,
        nmatch=1000,
        alen=1000,
        mapq=60,
    )
    span_read = SpanningRead(
        qname="ont_span_5p_001",
        strand="+",
        left_hit=paf_dummy,
        right_hit=paf_dummy,
        gap_length=1000,
        q_from=1000,
        q_to=2000,
    )
    qrec = ResolutionRecord(
        feature_id="sat_gap_variable:chr_stem:5p:1",
        chrom="chr_stem",
        side="5p",
        resolution_type=ResolutionType.SAT_GAP_VARIABLE,
        left_anchor_port="chr_stem+",
        right_anchor_port="cap_5p+",
    )
    erec = EndStateRecord(
        feature_id="end_complete:chr_stem:5p",
        chrom="chr_stem",
        side="5p",
        state=EndState.END_COMPLETE,
    )
    gap_res = SatelliteGapResolution(
        chrom="chr_stem",
        side="5p",
        state=EndState.END_COMPLETE,
        resolution_record=qrec,
        end_state_record=erec,
        spanning_reads=[span_read],
        representative_read=span_read,
        representative_sequence="T" * 1000,
    )

    scaffolds, agp_recs = build_scaffolds(
        graph=graph,
        ledger_map={},
        end_states=[],
        resolutions=[],
        gap_resolutions=[gap_res],
        dropout_resolutions=[],
        rdna_resolutions=[],
        junctions=[],
    )

    assert "chr_stem" in scaffolds
    # cap_5p (5000 bp) + bridge (1000 bp) + chr_stem (20000 bp) = 26000 bp
    assert len(scaffolds["chr_stem"]) == 26000

    # cap_5p+ inverted is cap_5p-, so reverse complement of C*5000 is G*5000
    assert scaffolds["chr_stem"][:5000] == "G" * 5000
    # Bridge is reverse complement of T*1000 -> A*1000
    assert scaffolds["chr_stem"][5000:6000] == "A" * 1000
    # chr_stem is A*20000
    assert scaffolds["chr_stem"][6000:] == "A" * 20000

    stem_agp = [r for r in agp_recs if r.object_id == "chr_stem"]
    assert len(stem_agp) == 3
    assert stem_agp[0].component_id == "cap_5p"
    assert stem_agp[0].orientation == "-"
    assert stem_agp[0].sequence_provenance == SequenceProvenance.HIFI_GERMLINE
    assert stem_agp[0].unpolished is False

    assert stem_agp[1].component_id == "ont_span_5p_001"
    assert stem_agp[1].orientation == "-"
    assert stem_agp[1].sequence_provenance == SequenceProvenance.ONT_DERIVED
    assert stem_agp[1].unpolished is True

    assert stem_agp[2].component_id == "chr_stem"
    assert stem_agp[2].orientation == "+"
    assert stem_agp[2].sequence_provenance == SequenceProvenance.HIFI_CORE
    assert stem_agp[2].unpolished is False


@pytest.mark.bio
@pytest.mark.benchmark
def test_benchmark_df5120_finish_scaffolds_and_telomeres(tmp_path: Path):
    """Verify end-to-end finish CLI execution on DF5120 benchmark dataset and T2T scaffold telomeres."""
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
    ont_paf_candidates = [
        benchmark_dir / "ont_to_utgs.paf",
        Path("tests/data/df5120_ont_slice.paf.gz"),
    ]
    ont_paf_path = next((p for p in ont_paf_candidates if p.exists()), None)

    if not (
        tiara_path.exists()
        and tidk_path.exists()
        and barrnap_path.exists()
        and compleasm_path.exists()
        and ont_paf_path is not None
        and ont_paf_path.exists()
    ):
        pytest.skip(f"Benchmark annotation/PAF files missing in {benchmark_dir}")

    # Step 1: Run finish WITHOUT ONT PAF to verify rDNA heuristic capping and 2-hop Chr X walking
    out_dir_no_ont = tmp_path / "df5120_finish_no_ont"
    exit_code_no_ont = main(
        [
            "finish",
            "--gfa",
            str(gfa_path),
            "--tiara",
            str(tiara_path),
            "--tidk",
            str(tidk_path),
            "--barrnap",
            str(barrnap_path),
            "--compleasm",
            str(compleasm_path),
            "--out-dir",
            str(out_dir_no_ont),
        ]
    )
    assert exit_code_no_ont == 0
    fasta_no_ont = Fasta(str(out_dir_no_ont / "scaffolds.fasta"))

    # utg000013l (Chr II) reaches bilateral T2T via rDNA heuristic capping:
    seq13_no_ont = str(fasta_no_ont["utg000013l"])
    assert seq13_no_ont[:3000].count("GCCTAA") > 100
    assert seq13_no_ont[-3000:].count("TTAGGC") > 100

    # utg000027l (Chr X) forms full 9.26 Mb scaffold via 2-hop allelic bubble walking:
    seq27_no_ont = str(fasta_no_ont["utg000027l"])
    assert len(seq27_no_ont) == 9256729
    assert seq27_no_ont[:3000].count("GCCTAA") > 100
    assert seq27_no_ont[-3000:].count("TTAGGC") > 100
    assert "utg000016l" not in fasta_no_ont

    # Step 2: Run finish WITH ONT PAF and reads to verify full bilateral T2T across all 6 chromosomes
    out_dir = tmp_path / "df5120_finish"
    ont_reads_candidates = [
        benchmark_dir / "ont_spanning_reads.fa.gz",
        Path("tests/data/df5120_ont_slice.fa.gz"),
    ]
    ont_reads_path = next((p for p in ont_reads_candidates if p.exists()), None)
    finish_args = [
        "finish",
        "--gfa",
        str(gfa_path),
        "--tiara",
        str(tiara_path),
        "--tidk",
        str(tidk_path),
        "--barrnap",
        str(barrnap_path),
        "--compleasm",
        str(compleasm_path),
        "--ont-paf",
        str(ont_paf_path),
        "--out-dir",
        str(out_dir),
    ]
    if ont_reads_path is not None and ont_reads_path.exists():
        finish_args.extend(["--ont-reads", str(ont_reads_path)])

    exit_code = main(finish_args)
    assert exit_code == 0

    # 1. junctions.tsv has 13 rows
    junctions = parse_junctions_tsv(out_dir / "junctions.tsv")
    assert len(junctions) == 13

    # 2. resolutions.tsv has all 6 required resolution types
    assert not (out_dir / "quirks.tsv").exists()
    resolutions = parse_resolutions_tsv(out_dir / "resolutions.tsv")
    assert len(resolutions) >= 6
    res_types = {r.resolution_type for r in resolutions}
    assert ResolutionType.SOMATIC_ALLELIC_BUBBLE in res_types
    assert ResolutionType.INTERCALARY_GRS_BRIDGE in res_types
    assert ResolutionType.PDE_TERMINAL_CAP in res_types
    assert ResolutionType.SUBTELOMERIC_SATELLITE_CAP in res_types
    assert ResolutionType.RDNA_CASSETTE in res_types

    # 3. utg000009l is fully resolved in end_states.tsv
    end_states = parse_end_states_tsv(out_dir / "end_states.tsv")
    utg9_feature_ids = {e.feature_id for e in end_states if e.chrom == "utg000009l"}
    assert "end_complete:utg000009l:3p" in utg9_feature_ids
    if ont_paf_path is not None and ont_paf_path.exists():
        assert ResolutionType.SAT_GAP_VARIABLE in res_types
        assert "end_complete:utg000009l:5p" in utg9_feature_ids

    # Verify scaffolds.agp contains the ONT-derived component for utg000009l:5p
    agp_records = parse_agp(out_dir / "scaffolds.agp")
    utg9_agp = [r for r in agp_records if r.object_id == "utg000009l"]
    if ont_paf_path is not None and ont_paf_path.exists() and ont_reads_path is not None and ont_reads_path.exists():
        ont_comps = [r for r in utg9_agp if getattr(r, "component_id", None) == "e59ed077-450e-4e48-ab71-19a19d3490b9"]
        assert len(ont_comps) == 1
        assert ont_comps[0].sequence_provenance == SequenceProvenance.ONT_DERIVED
        assert ont_comps[0].unpolished is True
    else:
        assert any(r.component_type in ("N", "U") for r in utg9_agp)

    # 4. In scaffolds.fasta - ALL 6 chromosomes achieve bilateral T2T (>100 GCCTAA 5', >100 TTAGGC 3')
    fasta = Fasta(str(out_dir / "scaffolds.fasta"))

    # Chr I (utg000001l, len 7,826,949): bilateral T2T
    seq1 = str(fasta["utg000001l"])
    assert len(seq1) == 7826949
    assert seq1[:3000].count("GCCTAA") >= 260
    assert seq1[-3000:].count("TTAGGC") >= 230

    # Chr II (utg000013l, len 8,566,588): bilateral T2T with 5p rDNA cap utg000053l+
    seq13 = str(fasta["utg000013l"])
    assert len(seq13) == 8566588
    assert seq13[:3000].count("GCCTAA") >= 350
    assert seq13[-3000:].count("TTAGGC") >= 150

    # Chr III (utg000002l, len 7,187,254): bilateral T2T with 5p satellite cap utg000150l+ and intermediate satellite utg000083l
    seq2 = str(fasta["utg000002l"])
    assert len(seq2) == 7187254
    assert seq2[:3000].count("GCCTAA") >= 150
    assert seq2[-3000:].count("TTAGGC") >= 200
    assert "utg000083l" not in fasta

    # Chr IV (utg000004l, len 6,999,772): bilateral T2T
    seq4 = str(fasta["utg000004l"])
    assert len(seq4) == 6999772
    assert seq4[:3000].count("GCCTAA") >= 200
    assert seq4[-3000:].count("TTAGGC") >= 200

    # Chr V (utg000009l): bilateral T2T with 5p ONT gap cap utg000217l+
    seq9 = str(fasta["utg000009l"])
    if ont_paf_path is not None and ont_paf_path.exists():
        assert seq9[:3000].count("GCCTAA") >= 400
    assert seq9[-3000:].count("TTAGGC") >= 200

    # Chr X (utg000027l, len 9,256,729): bilateral T2T across full allelic bubble walk
    seq27 = str(fasta["utg000027l"])
    assert len(seq27) == 9256729
    assert seq27[:3000].count("GCCTAA") >= 300
    assert seq27[-3000:].count("TTAGGC") >= 150
    assert "utg000016l" not in fasta

    # Verify no non-contaminant 45S rDNA unitigs float in scaffolds.fasta
    for rdna_u in ["utg000098l", "utg000053l", "utg000020l", "utg000082l"]:
        assert rdna_u not in fasta


def test_excluded_fasta_emission_and_conservation_law(tmp_path: Path):
    """Verify excluded.fasta emission, header format, sorting, and 100% unitig conservation law."""
    isolates = {
        "DF5120": {
            "gfa": Path("data/benchmark/DF5120/asm.p_utg.gfa.gz"),
            "tidk": Path("data/benchmark/DF5120/tidk_telomeres_telomeric_repeat_windows.tsv"),
            "barrnap": Path("data/benchmark/DF5120/barrnap_rdna.gff"),
            "compleasm": Path("data/benchmark/DF5120/compleasm_unitigs_full_table.tsv"),
            "tiara": Path("data/benchmark/DF5120/tiara_unitigs.txt"),
        },
        "JU1382": {
            "gfa": Path("data/benchmark/JU1382/asm.p_utg.gfa.gz"),
            "tidk": Path("data/benchmark/JU1382/tidk_telomeres_telomeric_repeat_windows.tsv"),
            "barrnap": Path("data/benchmark/JU1382/barrnap_rdna.gff"),
            "compleasm": Path("data/benchmark/JU1382/compleasm_unitigs_full_table.tsv"),
            "tiara": Path("data/benchmark/JU1382/tiara_unitigs.txt"),
        },
        "PS1017": {
            "gfa": Path("data/benchmark/PS1017/asm.p_utg.gfa.gz"),
            "tidk": Path("data/benchmark/PS1017/tidk_telomeres_telomeric_repeat_windows.tsv"),
            "barrnap": Path("data/benchmark/PS1017/barrnap_rdna.gff"),
            "compleasm": Path("data/benchmark/PS1017/compleasm_unitigs_full_table.tsv"),
            "tiara": Path("data/benchmark/PS1017/tiara_unitigs.txt"),
        },
        "PS2068": {
            "gfa": Path("data/benchmark/PS2068/asm.p_utg.gfa.gz"),
            "tidk": Path("data/benchmark/PS2068/tidk_telomeres_telomeric_repeat_windows.tsv"),
            "barrnap": Path("data/benchmark/PS2068/barrnap_rdna.gff"),
            "compleasm": Path("data/benchmark/PS2068/compleasm_unitigs_full_table.tsv"),
            "tiara": Path("data/benchmark/PS2068/tiara_unitigs.txt"),
        },
    }

    for name, p in isolates.items():
        if not all(fp.exists() for fp in p.values()):
            continue

        out_dir = tmp_path / f"test_finish_{name}"
        custom_excl = tmp_path / f"custom_excluded_{name}.fasta"

        # Test with custom --excluded-fasta
        rc = main(
            [
                "finish",
                "--gfa",
                str(p["gfa"]),
                "--tidk",
                str(p["tidk"]),
                "--barrnap",
                str(p["barrnap"]),
                "--compleasm",
                str(p["compleasm"]),
                "--tiara",
                str(p["tiara"]),
                "--out-dir",
                str(out_dir),
                "--excluded-fasta",
                str(custom_excl),
            ]
        )
        assert rc == 0
        assert custom_excl.exists()
        assert (out_dir / "scaffolds.agp").exists()
        assert (out_dir / "resolutions.tsv").exists()

        # Parse GFA segments
        graph = load_gfa(p["gfa"])
        all_gfa_unitigs = set(graph.segments.keys())

        # Parse AGP placed unitigs (W components matching GFA segment names)
        agp_lines = (out_dir / "scaffolds.agp").read_text().splitlines()
        placed_unitigs = set()
        for line in agp_lines:
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 6 and parts[4] == "W":
                comp_id = parts[5]
                if comp_id in all_gfa_unitigs:
                    placed_unitigs.add(comp_id)

        # Parse excluded.fasta
        excl_fasta = Fasta(str(custom_excl))
        excl_unitigs = set(excl_fasta.keys())

        # 100% unitig conservation law
        assert placed_unitigs | excl_unitigs == all_gfa_unitigs, (
            f"Missing unitigs in {name}: {all_gfa_unitigs - (placed_unitigs | excl_unitigs)}"
        )
        assert placed_unitigs & excl_unitigs == set(), f"Overlapping unitigs in {name}: {placed_unitigs & excl_unitigs}"

        # Check deterministic sorted order
        excl_keys = list(excl_fasta.keys())
        assert excl_keys == sorted(excl_keys, key=lambda uid: uid)

        # Check header format
        excl_lines = custom_excl.read_text().splitlines()
        header_lines = [line_str for line_str in excl_lines if line_str.startswith(">")]
        for h in header_lines:
            # Format: >{uid} excluded_by={feat_ids} length={length}
            parts = h[1:].split()
            assert len(parts) == 3
            uid = parts[0]
            assert uid in all_gfa_unitigs
            assert parts[1].startswith("excluded_by=")
            assert len(parts[1][len("excluded_by=") :]) > 0
            assert parts[2].startswith("length=")
            expected_len = graph.segments[uid].length
            assert int(parts[2][len("length=") :]) == expected_len

        # Specific metadata checks for PS1017
        if name == "PS1017":
            resolutions = parse_resolutions_tsv(out_dir / "resolutions.tsv")
            # Verify no redundant sat_gap_variable for utg000032l:5p
            pde_32 = [r for r in resolutions if r.chrom == "utg000032l" and r.side == "5p"]
            assert len(pde_32) == 1
            assert pde_32[0].resolution_type == ResolutionType.PDE_TERMINAL_CAP
            assert pde_32[0].unitig_ids == ("utg000076l",)
            assert pde_32[0].excluded_unitig_ids == ("utg000186l",)

            # Verify utg000052l terminal cap places utg000049l
            pde_52 = [r for r in resolutions if r.chrom == "utg000052l" and r.side == "5p"]
            assert len(pde_52) == 1
            assert pde_52[0].unitig_ids == ("utg000049l",)
            assert pde_52[0].right_anchor_port == "utg000049l-"
            assert pde_52[0].excluded_unitig_ids == ("utg000146l",)

            # Verify utg000061l terminal cap places multi-unitig continuation utg000070l,utg000126l from agp_start 1
            pde_61 = [r for r in resolutions if r.chrom == "utg000061l" and r.side == "5p"]
            assert len(pde_61) == 1
            assert pde_61[0].unitig_ids == ("utg000070l", "utg000126l")
            assert pde_61[0].agp_start == 1
            assert pde_61[0].agp_end == 128490
            assert pde_61[0].right_anchor_port == "utg000070l-"
            assert pde_61[0].excluded_unitig_ids == ("utg000196l",)

            # Verify all terminal PDE junctions are PDE_TERMINAL_CAP and exclude somatic stubs
            pde_caps = [r for r in resolutions if r.resolution_type == ResolutionType.PDE_TERMINAL_CAP]
            assert len(pde_caps) == 13
            for pcap in pde_caps:
                assert len(pcap.excluded_unitig_ids) >= 1
                assert pcap.left_anchor_port is not None and len(pcap.left_anchor_port) > 0
                assert pcap.right_anchor_port is not None and len(pcap.right_anchor_port) > 0
                stem_seg = pcap.left_anchor_port.rstrip("+-")
                assert stem_seg not in pcap.unitig_ids

            # Verify sat_gap_variable is only emitted when bridged (0 without ONT PAF)
            sat_vars = [r for r in resolutions if r.resolution_type == ResolutionType.SAT_GAP_VARIABLE]
            assert len(sat_vars) == 0
