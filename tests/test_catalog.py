from __future__ import annotations

import csv
import hashlib
import json
import random
from pathlib import Path
from types import SimpleNamespace

import pytest

from gfapdefinisher.catalog import (
    EndCategory,
    HorizonReason,
    build_end_catalog,
    emit_end_catalog,
    load_end_evidence,
)
from gfapdefinisher.cli import main
from gfapdefinisher.emitter import AgpGapRow, AgpRecord
from gfapdefinisher.graph import BidirectedGraph


def _random_seq(rng: random.Random, n: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(n))


def test_ends_sharing_an_unresolved_satellite_form_one_group_with_a_fusion_interpretation():
    rng = random.Random(7)
    unit = _random_seq(rng, 132)
    while "TTAGGC" in unit * 2 or "GCCTAA" in unit * 2:
        unit = _random_seq(rng, 132)

    graph = BidirectedGraph()
    seqs = {
        "coreA": _random_seq(rng, 20000) + unit * 20,
        "coreB": _random_seq(rng, 18000) + unit * 20,
        "capL": _random_seq(rng, 3000),
        "capM": _random_seq(rng, 3000),
        "capX": unit * 10 + _random_seq(rng, 2000),
        "capY": unit * 8 + _random_seq(rng, 2000),
        "capW": _random_seq(rng, 4000),
    }
    for name, seq in seqs.items():
        graph.add_segment(name, len(seq), depth=60.0 if name.startswith("core") else 5.0, sequence=seq)
    graph.add_edge("capL+", "coreA+", "100M")
    graph.add_edge("capM+", "coreB+", "100M")

    telomeric = {"capL-", "capM-", "capX+", "capY+", "capW+"}
    ledger = {
        f"{name}{s}": SimpleNamespace(has_terminal_telomere=f"{name}{s}" in telomeric, rdna_hit_count=0)
        for name in seqs
        for s in "+-"
    }
    orphan_caps = {name: (f"{name}+", f"{name}-") for name in ("capX", "capY", "capW")}

    agp: list[AgpRecord | AgpGapRow] = []
    for obj, (cap5, core, cap3) in {"SUPER_1": ("capL", "coreA", "capX"), "SUPER_2": ("capM", "coreB", "capY")}.items():
        pos = 1
        for part, (kind, name) in enumerate([("W", cap5), ("W", core), ("U", None), ("W", cap3)], start=1):
            if kind == "U":
                agp.append(AgpGapRow(obj, pos, pos + 99, part, "U", 100, "repeat", "yes", "satellite-kmer"))
                pos += 100
                continue
            length = len(seqs[name])
            agp.append(AgpRecord(obj, pos, pos + length - 1, part, "W", name, 1, length, "+"))
            pos += length

    entries, groups = build_end_catalog(graph, ledger, agp, orphan_caps, core_coverage=60.0)
    by_end = {e.end_id: e for e in entries}

    assert by_end["SUPER_1:5p"].category == EndCategory.RESOLVED
    assert by_end["SUPER_2:5p"].category == EndCategory.RESOLVED
    for end_id, core, chosen in (("SUPER_1:3p", "coreA+", "capX+"), ("SUPER_2:3p", "coreB+", "capY+")):
        entry = by_end[end_id]
        assert entry.horizon_port == core
        assert entry.horizon_reason == HorizonReason.SATELLITE
        assert entry.tip_repeat_period_bp == 132
        assert entry.finish_choice == chosen
        assert {c.port for c in entry.candidates} == {"capX+", "capY+"}
        assert entry.category == EndCategory.SHARED

    shared = next(g for g in groups if g.ends)
    assert shared.ends == ("SUPER_1:3p", "SUPER_2:3p")
    assert shared.fusion_interpretations == ("SUPER_1:3p<->SUPER_2:3p",)
    assert set(shared.shared_unitigs) == {"capX", "capY"}
    assert [g.candidate_endings for g in groups if not g.ends] == [("capW+",)]


def _make_end_evidence_payload(primary_sha: str) -> dict:
    return {
        "version": 1,
        "coordinate_system": "zero-based half-open on forward sequences",
        "graph_sha256": {"primary": primary_sha, "raw": "deadbeef"},
        "families": [
            {
                "cluster_id": "fam_A",
                "members": [
                    {
                        "end_id": "mol_1:3p",
                        "molecule_id": "mol_1",
                        "read_end": "3p",
                        "tract_start": 10,
                        "tract_end": 210,
                    },
                    {
                        "end_id": "mol_1:5p",
                        "molecule_id": "mol_1",
                        "read_end": "5p",
                        "tract_start": 0,
                        "tract_end": 180,
                    },
                    {
                        "end_id": "mol_2:3p",
                        "molecule_id": "mol_2",
                        "read_end": "3p",
                        "tract_start": 5,
                        "tract_end": 195,
                    },
                ],
                "haplotype_patterns": [{"pattern_id": "pat_1", "read_count": 2, "alleles": "ACGT"}],
            },
            {
                "cluster_id": "fam_raw_only",
                "members": [
                    {
                        "end_id": "mol_raw:3p",
                        "molecule_id": "mol_raw",
                        "read_end": "3p",
                        "tract_start": 0,
                        "tract_end": 100,
                    }
                ],
                "haplotype_patterns": [],
            },
        ],
        "consensus_placements": [
            {
                "cluster_id": "fam_A",
                "source": "primary",
                "record": "coreA",
                "strand": "+",
                "start": 18000,
                "end": 18500,
                "target_bp": 22000,
                "consensus_start": 0,
                "consensus_end": 500,
                "matches": 495,
                "alignment_bp": 500,
                "mapq": 60,
                "family_molecules": 2,
                "patterns": "pat_1",
            },
            {
                "cluster_id": "fam_A",
                "source": "primary",
                "record": "capX",
                "strand": "+",
                "start": 0,
                "end": 500,
                "target_bp": 3000,
                "consensus_start": 0,
                "consensus_end": 500,
                "matches": 490,
                "alignment_bp": 500,
                "mapq": 50,
                "family_molecules": 2,
                "patterns": "pat_1",
            },
            {
                "cluster_id": "fam_A",
                "source": "primary",
                "record": "capY",
                "strand": "+",
                "start": 0,
                "end": 500,
                "target_bp": 3000,
                "consensus_start": 0,
                "consensus_end": 500,
                "matches": 485,
                "alignment_bp": 500,
                "mapq": 45,
                "family_molecules": 2,
                "patterns": "pat_1",
            },
            {
                "cluster_id": "fam_raw_only",
                "source": "raw",
                "record": "capX",
                "strand": "+",
                "start": 0,
                "end": 400,
                "target_bp": 3000,
                "consensus_start": 0,
                "consensus_end": 400,
                "matches": 395,
                "alignment_bp": 400,
                "mapq": 60,
                "family_molecules": 1,
                "patterns": "",
            },
        ],
    }


def test_load_end_evidence_validation(tmp_path: Path):
    gfa = tmp_path / "primary.gfa"
    gfa.write_bytes(b"H\tVN:Z:1.0\n")
    primary_sha = hashlib.sha256(gfa.read_bytes()).hexdigest()

    ev_path = tmp_path / "end_evidence.json"
    payload = _make_end_evidence_payload(primary_sha)
    ev_path.write_text(json.dumps(payload), encoding="utf-8")

    loaded = load_end_evidence(ev_path, primary_gfa=gfa)
    assert loaded.version == 1
    assert len(loaded.families) == 2
    assert len(loaded.consensus_placements) == 4

    bad_ver = {**payload, "version": 2}
    ev_path.write_text(json.dumps(bad_ver), encoding="utf-8")
    with pytest.raises(ValueError, match="Unsupported end-evidence version"):
        load_end_evidence(ev_path, primary_gfa=gfa)

    bad_coord = {**payload, "coordinate_system": "one-based closed"}
    ev_path.write_text(json.dumps(bad_coord), encoding="utf-8")
    with pytest.raises(ValueError, match="Unsupported end-evidence coordinate_system"):
        load_end_evidence(ev_path, primary_gfa=gfa)

    bad_sha = {**payload, "graph_sha256": {"primary": "0" * 64}}
    ev_path.write_text(json.dumps(bad_sha), encoding="utf-8")
    with pytest.raises(ValueError, match="Primary GFA SHA-256 mismatch"):
        load_end_evidence(ev_path, primary_gfa=gfa)


def test_end_catalog_attaches_families_and_molecules(tmp_path: Path):
    rng = random.Random(11)
    graph = BidirectedGraph()
    seqs = {
        "capL": _random_seq(rng, 3000),
        "coreA": _random_seq(rng, 20000),
        "capX": _random_seq(rng, 3000),
        "capY": _random_seq(rng, 3000),
    }
    for name, seq in seqs.items():
        graph.add_segment(name, len(seq), depth=60.0 if name == "coreA" else 5.0, sequence=seq)
    graph.add_edge("capL+", "coreA+", "100M")

    telomeric = {"capL-", "capX+", "capY+"}
    ledger = {
        f"{name}{s}": SimpleNamespace(has_terminal_telomere=f"{name}{s}" in telomeric, rdna_hit_count=0)
        for name in seqs
        for s in "+-"
    }
    orphan_caps = {name: (f"{name}+", f"{name}-") for name in ("capX", "capY")}
    agp = [
        AgpRecord("SUPER_1", 1, len(seqs["capL"]), 1, "W", "capL", 1, len(seqs["capL"]), "+"),
        AgpRecord(
            "SUPER_1",
            len(seqs["capL"]) + 1,
            len(seqs["capL"]) + len(seqs["coreA"]),
            2,
            "W",
            "coreA",
            1,
            len(seqs["coreA"]),
            "+",
        ),
    ]

    gfa = tmp_path / "primary.gfa"
    gfa.write_bytes(b"H\tVN:Z:1.0\n")
    ev_path = tmp_path / "end_evidence.json"
    ev_path.write_text(json.dumps(_make_end_evidence_payload(hashlib.sha256(gfa.read_bytes()).hexdigest())))
    end_evidence = load_end_evidence(ev_path, primary_gfa=gfa)

    entries, groups = build_end_catalog(
        graph,
        ledger,
        agp,
        orphan_caps,
        core_coverage=60.0,
        end_evidence=end_evidence,
    )
    by_end = {e.end_id: e for e in entries}
    end_3p = by_end["SUPER_1:3p"]
    # Both competing primary placements (capX+ and capY+) are retained; raw-only fam_raw_only is ignored.
    assert [c.port for c in end_3p.candidates] == ["capX+", "capY+"]
    for cand in end_3p.candidates:
        assert cand.evidence == ("family:fam_A",)
        assert cand.family_ids == ("fam_A",)
        assert cand.molecule_ids == ("mol_1", "mol_2")

    emit_end_catalog(entries, groups, tmp_path)
    with (tmp_path / "end_catalog.tsv").open() as fh:
        rows = {r["end_id"]: r for r in csv.DictReader(fh, delimiter="\t")}
    assert rows["SUPER_1:3p"]["candidate_families"] == "fam_A;fam_A"
    assert rows["SUPER_1:3p"]["candidate_molecules"] == "mol_1,mol_2;mol_1,mol_2"
    assert rows["SUPER_1:5p"]["candidate_families"] == ""
    assert rows["SUPER_1:5p"]["candidate_molecules"] == ""


def test_cli_finish_end_evidence_flag(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    with pytest.raises(SystemExit):
        main(["finish", "--help"])
    assert "--end-evidence" in capsys.readouterr().out

    gfa = tmp_path / "asm.gfa"
    core = "A" * 120000 + "TTAGGC" * 35
    cap_5p = "GCCTAA" * 35 + "A" * 10000
    som_5p = "GCCTAA" * 35 + "A" * 10000
    cap_orphan = "G" * 10000 + "TTAGGC" * 35
    gfa.write_text(
        f"H\tVN:Z:1.0\n"
        f"S\tcap5\t{cap_5p}\tLN:i:{len(cap_5p)}\tdp:f:6.0\n"
        f"S\tsom5\t{som_5p}\tLN:i:{len(som_5p)}\tdp:f:30.0\n"
        f"S\tcore1\t{core}\tLN:i:{len(core)}\tdp:f:60.0\n"
        f"S\tcap3\t{cap_orphan}\tLN:i:{len(cap_orphan)}\tdp:f:5.0\n"
        f"L\tcore1\t-\tcap5\t-\t100M\n"
        f"L\tcore1\t-\tsom5\t-\t100M\n",
        encoding="utf-8",
    )
    primary_sha = hashlib.sha256(gfa.read_bytes()).hexdigest()
    ev_path = tmp_path / "phased_end_evidence.json"
    ev_path.write_text(
        json.dumps(
            {
                "version": 1,
                "coordinate_system": "zero-based half-open on forward sequences",
                "graph_sha256": {"primary": primary_sha},
                "families": [
                    {
                        "cluster_id": "fam_7",
                        "members": [
                            {
                                "end_id": "read_x:3p",
                                "molecule_id": "read_x",
                                "read_end": "3p",
                                "tract_start": 0,
                                "tract_end": 120,
                            }
                        ],
                        "haplotype_patterns": [],
                    }
                ],
                "consensus_placements": [
                    {
                        "cluster_id": "fam_7",
                        "source": "primary",
                        "record": "core1",
                        "strand": "+",
                        "start": len(core) - 200,
                        "end": len(core),
                        "target_bp": len(core),
                        "consensus_start": 0,
                        "consensus_end": 200,
                        "matches": 198,
                        "alignment_bp": 200,
                        "mapq": 60,
                        "family_molecules": 1,
                        "patterns": "",
                    },
                    {
                        "cluster_id": "fam_7",
                        "source": "primary",
                        "record": "cap3",
                        "strand": "+",
                        "start": len(cap_orphan) - 200,
                        "end": len(cap_orphan),
                        "target_bp": len(cap_orphan),
                        "consensus_start": 0,
                        "consensus_end": 200,
                        "matches": 198,
                        "alignment_bp": 200,
                        "mapq": 60,
                        "family_molecules": 1,
                        "patterns": "",
                    },
                ],
            }
        ),
        encoding="utf-8",
    )

    out_dir = tmp_path / "out"
    assert (
        main(
            [
                "finish",
                "--gfa",
                str(gfa),
                "--out-dir",
                str(out_dir),
                "--nominate-super",
                "--end-evidence",
                str(ev_path),
            ]
        )
        == 0
    )
    with (out_dir / "end_catalog.tsv").open() as fh:
        rows = {r["end_id"]: r for r in csv.DictReader(fh, delimiter="\t")}
    assert rows["SUPER_1:3p"]["candidate_endings"] == "cap3+[family:fam_7]"
    assert rows["SUPER_1:3p"]["candidate_families"] == "fam_7"
    assert rows["SUPER_1:3p"]["candidate_molecules"] == "read_x"
