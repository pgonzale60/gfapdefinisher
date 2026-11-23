from __future__ import annotations

import random
from types import SimpleNamespace

from gfapdefinisher.catalog import EndCategory, HorizonReason, build_end_catalog
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
