"""Tests for the gfapdefinisher command-line interface."""

from __future__ import annotations

from pathlib import Path

import pytest

from gfapdefinisher.cli import main
from gfapdefinisher.policy import DEFAULT_POLICY
from gfapdefinisher.quirks import (
    EndState,
    IntervalQuirk,
    parse_end_states_tsv,
    parse_junctions_tsv,
    parse_quirks_tsv,
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
    quirks_file = out_dir / "quirks.tsv"

    assert ledger_file.is_file()
    assert junctions_file.is_file()
    assert end_states_file.is_file()
    assert quirks_file.is_file()

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

    quirk_recs = parse_quirks_tsv(quirks_file)
    assert len(quirk_recs) >= 1
    assert any(q.quirk_type == IntervalQuirk.SAT_GAP_VARIABLE for q in quirk_recs)


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
    quirks_file = out_dir / "quirks.tsv"

    end_state_recs = parse_end_states_tsv(end_states_file)
    # The gap should now be resolved to END_COMPLETE
    complete_ends = [es for es in end_state_recs if es.state == EndState.END_COMPLETE]
    assert len(complete_ends) >= 1
    assert any(es.resolved is True for es in complete_ends)
    resolved_end = [es for es in complete_ends if es.chrom == "chr_stem"][0]
    assert resolved_end.somatic_port == "chr_som+"
    assert resolved_end.germline_port == "orphan_cap+"
    assert not any(es.chrom == "orphan_cap" for es in end_state_recs)

    quirk_recs = parse_quirks_tsv(quirks_file)
    sat_quirks = [q for q in quirk_recs if q.quirk_type == IntervalQuirk.SAT_GAP_VARIABLE]
    assert len(sat_quirks) >= 1
    bridged_quirk = sat_quirks[0]
    assert bridged_quirk.spanning_read_count == 3
    assert bridged_quirk.representative_sequence_selected is True
    assert bridged_quirk.agp_linkage_evidence == "ont-spanning"
    assert bridged_quirk.length_median_bp == 1500


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
    quirks_file = out_dir / "quirks.tsv"
    quirk_recs = parse_quirks_tsv(quirks_file)
    sat_quirks = [q for q in quirk_recs if q.quirk_type == IntervalQuirk.SAT_GAP_VARIABLE]
    assert len(sat_quirks) >= 1
    assert sat_quirks[0].representative_sequence_selected is True


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
    assert ">chr_som" in fasta_content

    agp_content = agp_file.read_text(encoding="utf-8")
    assert agp_content.startswith("##agp-version 2.1\n")
    assert "chr_stem" in agp_content
    assert "chr_som" in agp_content


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
    assert ">TEST_chr_som" in fa_text

    agp_text = custom_agp.read_text(encoding="utf-8")
    assert "TEST_chr_stem" in agp_text


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

    quirks = parse_quirks_tsv(out_dir / "quirks.tsv")
    rdna_quirks = [q for q in quirks if q.quirk_type == IntervalQuirk.RDNA_CASSETTE]
    assert len(rdna_quirks) >= 1
    assert rdna_quirks[0].chrom == "chr_arm"
    assert "distal cap distal_cap+" in rdna_quirks[0].notes

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
