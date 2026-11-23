"""Regression coverage for replaying reviewed joins and the DF5120 milestone."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from gfapdefinisher.cli import main
from gfapdefinisher.curation import curate

ROOT = Path(__file__).resolve().parent.parent
RECIPE = ROOT / "benchmarks" / "df5120_t2t.json"


def test_curation_checks_read_anchors_and_keeps_a_port_gap(tmp_path: Path):
    left = "GCCTAA" * 32 + "AACG"
    right = "AGTC" + "TTAGGC" * 32
    graph = tmp_path / "source.gfa"
    graph.write_text(f"S\tleft\t{left}\nS\tright\t{right}\n")
    ont = tmp_path / "reads.fa"
    ont.write_text(">bridge\nAACTTCTAGTC\n")
    path = tmp_path / "recipe.json"
    left_part = {"source": "primary", "record": "left"}
    right_part = {"source": "primary", "record": "right"}
    chromosome = {"cap_5p": left_part, "cap_3p": right_part}
    path.write_text(
        json.dumps(
            {
                "version": 1,
                "gap_bp": 200,
                "bridge_anchor_bp": 4,
                "cap_match_bp": 20,
                "line_width": 60,
                "chromosomes": [
                    {
                        **chromosome,
                        "name": "bridged",
                        "gap_count": 0,
                        "parts": [
                            left_part,
                            {"source": "ont", "record": "bridge", "start": 4, "end": 7},
                            right_part,
                        ],
                    },
                    {
                        **chromosome,
                        "name": "gapped",
                        "gap_count": 1,
                        "parts": [left_part, {"gap": "satellite_port"}, right_part],
                    },
                ],
            }
        )
    )
    output = tmp_path / "curated.fa"
    sources = {"primary": graph, "ont": ont}
    with pytest.raises(ValueError, match="does not match both splice anchors"):
        curate(path, sources, output)
    assert not output.exists()

    ont.write_text(">bridge\nAACGTCTAGTC\n")
    digest = curate(path, sources, output)
    records = {
        lines[0]: "".join(lines[1:]) for lines in (part.splitlines() for part in output.read_text().split(">")[1:])
    }
    assert records == {"bridged": left + "TCT" + right, "gapped": left + "N" * 200 + right}
    assert curate(path, sources, output) == digest


@pytest.mark.benchmark
def test_df5120_curation_reproduces_milestone(tmp_path: Path):
    data = ROOT / "data" / "benchmark" / "DF5120"
    files = {
        "--gfa": data / "asm.p_utg.gfa.gz",
        "--raw-gfa": data / "asm.r_utg.gfa.gz",
        "--tiara": data / "tiara_unitigs.txt",
        "--barrnap": data / "barrnap_rdna.gff",
        "--ont-paf": data / "ont_to_utgs.paf",
        "--ont-reads": data / "ont_spanning_reads.fa.gz",
    }
    missing = [path.name for path in files.values() if not path.is_file()]
    if missing:
        pytest.skip("Missing DF5120 milestone inputs: " + ", ".join(missing))

    finish_dir = tmp_path / "finish"
    finish = ["finish", "--nominate-super", "--out-dir", str(finish_dir)]
    for flag in ("--gfa", "--tiara", "--barrnap", "--ont-paf", "--ont-reads"):
        finish.extend((flag, str(files[flag])))
    assert main(finish) == 0

    output = tmp_path / "DF5120.fa"
    apply = [
        "curate",
        "--recipe",
        str(RECIPE),
        "--scaffolds",
        str(finish_dir / "scaffolds.fasta"),
        "--out",
        str(output),
    ]
    for flag in ("--gfa", "--raw-gfa", "--ont-reads"):
        apply.extend((flag, str(files[flag])))
    assert main(apply) == 0
    assert hashlib.sha256(output.read_bytes()).hexdigest() == json.loads(RECIPE.read_text())["expected_sha256"]
