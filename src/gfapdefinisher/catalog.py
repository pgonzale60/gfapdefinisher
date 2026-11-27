"""Chromosome-end uncertainty catalog.

For every chromosome end, separate what the graph and reads establish from what
`finish` chose heuristically:

- **Certain part**: components from the core outward while each join is a graph
  overlap or read-supplied sequence.
- **Horizon**: the last certain port. Beyond it the end enters sequence the data
  cannot place: an unresolved repeat gap, or a telomere on a core-depth unitig
  (the somatic end; any germline extension is absent from the graph).
- **Candidate endings**: orphan telomeric caps linked to the horizon by graph
  reachability, spanning reads, the tip satellite k-mer, or rDNA content. Caps that
  another end reaches through certain joins are resolved occurrences and link no end.
- **Ambiguity groups**: connected components of ends and the unitigs they share.
  Two or more ends in one group admit a fusion interpretation from sequence alone.
"""

from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter, deque
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from enum import StrEnum
from io import DEFAULT_BUFFER_SIZE
from pathlib import Path
from typing import Any

from gfapdefinisher.graph import BidirectedGraph
from gfapdefinisher.models import SegmentPort, Strand
from gfapdefinisher.policy import DEFAULT_POLICY, Policy
from gfapdefinisher.resolver import (
    SATELLITE_MIN_TIP_REPEATS,
    OntBridgingEngine,
    PafHit,
    _port_boundary_slack,
    reverse_complement,
    tip_satellite_kmer,
)

CHROMOSOME_PREFIX = "SUPER_"


class EndCategory(StrEnum):
    RESOLVED = "resolved"
    BOUNDED = "bounded"
    SHARED = "shared"
    OPEN = "open"


class HorizonReason(StrEnum):
    SATELLITE = "satellite"
    RDNA = "rdna"
    CORE_TELOMERE = "core_telomere"
    NO_TELOMERE = "no_telomere"
    NO_GRAPH_EDGE = "no_graph_edge"
    GAP = "gap"


_LINKAGE_REASON = {"satellite-kmer": HorizonReason.SATELLITE, "rna_tangle": HorizonReason.RDNA}
# Gap linkages whose order and orientation the data establish; only the sequence is absent.
_CERTAIN_GAP_JOIN = {"ont-spanning": "read-gap", "graph": "graph-gap", "graph-walk": "graph-gap"}


@dataclass(frozen=True, slots=True)
class CandidateEnding:
    """A telomeric cap the end may terminate in, with the evidence linking it."""

    port: str
    evidence: tuple[str, ...]
    family_ids: tuple[str, ...] = ()
    molecule_ids: tuple[str, ...] = ()

    def format(self) -> str:
        return f"{self.port}[{','.join(self.evidence)}]"


@dataclass(frozen=True, slots=True)
class EndEvidenceMember:
    """Observed read end belonging to a phased telomere family."""

    end_id: str
    molecule_id: str
    read_end: str
    tract_start: int
    tract_end: int


@dataclass(frozen=True, slots=True)
class HaplotypePattern:
    """Local aggregate haplotype pattern within a phased telomere family."""

    pattern_id: str
    read_count: int
    alleles: str


@dataclass(frozen=True, slots=True)
class EndEvidenceFamily:
    """Phased telomere family record keyed by cluster_id."""

    cluster_id: str
    members: tuple[EndEvidenceMember, ...]
    haplotype_patterns: tuple[HaplotypePattern, ...] = ()


@dataclass(frozen=True, slots=True)
class ConsensusPlacement:
    """Oriented placement of a family consensus on a primary or raw graph unitig."""

    cluster_id: str
    source: str
    record: str
    strand: str
    start: int
    end: int
    target_bp: int
    consensus_start: int
    consensus_end: int
    matches: int
    alignment_bp: int
    mapq: int
    family_molecules: int
    patterns: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class PhasedEndEvidence:
    """Validated version-1 phased telomere end-evidence payload."""

    version: int
    coordinate_system: str
    graph_sha256: Mapping[str, str]
    families: tuple[EndEvidenceFamily, ...]
    consensus_placements: tuple[ConsensusPlacement, ...]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(DEFAULT_BUFFER_SIZE), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_end_evidence(path: Path, primary_gfa: Path | None = None) -> PhasedEndEvidence:
    """Load and validate a version-1 phased end-evidence JSON file."""
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, Mapping):
        raise ValueError("End-evidence JSON must be an object")
    version = data.get("version")
    if isinstance(version, bool) or version != 1:
        raise ValueError(f"Unsupported end-evidence version: {version!r}")
    coord_sys = data.get("coordinate_system")
    if coord_sys != "zero-based half-open on forward sequences":
        raise ValueError(f"Unsupported end-evidence coordinate_system: {coord_sys!r}")
    graph_sha256 = data.get("graph_sha256")
    if not isinstance(graph_sha256, Mapping) or not isinstance(graph_sha256.get("primary"), str):
        raise ValueError("End-evidence JSON must include graph_sha256.primary")
    if primary_gfa is not None:
        actual_sha = _sha256(primary_gfa)
        if graph_sha256["primary"] != actual_sha:
            raise ValueError(
                f"Primary GFA SHA-256 mismatch: end-evidence expects {graph_sha256['primary']}, "
                f"but {primary_gfa} has {actual_sha}"
            )

    families: list[EndEvidenceFamily] = []
    for fam in data.get("families", []):
        members = tuple(
            EndEvidenceMember(
                end_id=str(m["end_id"]),
                molecule_id=str(m["molecule_id"]),
                read_end=str(m["read_end"]),
                tract_start=int(m["tract_start"]),
                tract_end=int(m["tract_end"]),
            )
            for m in fam.get("members", [])
        )
        patterns = tuple(
            HaplotypePattern(
                pattern_id=str(p["pattern_id"]),
                read_count=int(p["read_count"]),
                alleles=str(p["alleles"]),
            )
            for p in fam.get("haplotype_patterns", [])
        )
        families.append(
            EndEvidenceFamily(
                cluster_id=str(fam["cluster_id"]),
                members=members,
                haplotype_patterns=patterns,
            )
        )

    placements: list[ConsensusPlacement] = []
    for pl in data.get("consensus_placements", []):
        source = str(pl["source"])
        if source not in ("primary", "raw"):
            raise ValueError(f"Invalid consensus_placement source: {source!r}")
        strand = str(pl["strand"])
        if strand not in ("+", "-"):
            raise ValueError(f"Invalid consensus_placement strand: {strand!r}")
        raw_patterns = pl.get("patterns", "")
        pat_tuple = tuple(p.strip() for p in str(raw_patterns).split(",") if p.strip())
        placements.append(
            ConsensusPlacement(
                cluster_id=str(pl["cluster_id"]),
                source=source,
                record=str(pl["record"]),
                strand=strand,
                start=int(pl["start"]),
                end=int(pl["end"]),
                target_bp=int(pl["target_bp"]),
                consensus_start=int(pl["consensus_start"]),
                consensus_end=int(pl["consensus_end"]),
                matches=int(pl["matches"]),
                alignment_bp=int(pl["alignment_bp"]),
                mapq=int(pl["mapq"]),
                family_molecules=int(pl["family_molecules"]),
                patterns=pat_tuple,
            )
        )

    return PhasedEndEvidence(
        version=version,
        coordinate_system=coord_sys,
        graph_sha256=dict(graph_sha256),
        families=tuple(families),
        consensus_placements=tuple(placements),
    )


@dataclass(slots=True)
class EndCatalogEntry:
    """Certainty assessment for one chromosome end."""

    object_id: str
    side: str
    certain_ports: tuple[str, ...]
    certain_joins: tuple[str, ...]
    horizon_port: str
    horizon_reason: HorizonReason | None
    finish_choice: str | None = None
    finish_linkage: str | None = None
    tip_repeat_kmer: str | None = None
    tip_repeat_period_bp: int | None = None
    reads_at_horizon: int | None = None
    max_read_overhang_bp: int | None = None
    candidates: list[CandidateEnding] = field(default_factory=list)
    block_unitigs: tuple[str, ...] = ()
    group_id: str | None = None
    category: EndCategory = EndCategory.RESOLVED

    @property
    def end_id(self) -> str:
        return f"{self.object_id}:{self.side}"


@dataclass(slots=True)
class AmbiguityGroup:
    """Ends and candidate endings connected through shared unresolved unitigs."""

    group_id: str
    ends: tuple[str, ...]
    candidate_endings: tuple[str, ...]
    shared_unitigs: tuple[str, ...]

    @property
    def fusion_interpretations(self) -> tuple[str, ...]:
        return tuple(f"{a}<->{b}" for i, a in enumerate(self.ends) for b in self.ends[i + 1 :])


@dataclass(frozen=True, slots=True)
class _Item:
    kind: str  # "W" for a sequence component, otherwise the gap component type
    port: str | None = None
    linkage: str | None = None


def _port(sym: str) -> SegmentPort:
    return SegmentPort(sym[:-1], Strand.from_symbol(sym[-1]))


def _flip(sym: str) -> str:
    return sym[:-1] + ("-" if sym.endswith("+") else "+")


def _outward_items(records: Sequence[Any], side: str, graph: BidirectedGraph) -> list[_Item] | None:
    """Components from the object's longest unitig to the chromosome tip, in walk orientation.

    The longest unitig anchors the walk because depth alone cannot separate core
    unitigs from collapsed rDNA or satellite unitigs.
    """
    items = [
        _Item("W", port=f"{r.component_id}{r.orientation}")
        if r.component_type == "W"
        else _Item(r.component_type, linkage=str(r.orientation))
        for r in sorted(records, key=lambda r: r.part_number)
    ]
    if side == "5p":
        items = [_Item("W", port=_flip(it.port)) if it.port else it for it in reversed(items)]
    lengths = [
        (graph.segments[it.port[:-1]].length, i)
        for i, it in enumerate(items)
        if it.port and it.port[:-1] in graph.segments
    ]
    if not lengths:
        return None
    return items[max(lengths)[1] :]


def _tip_period(oriented_seq: str, kmer: str, window: int) -> int | None:
    tip = oriented_seq[-window:]
    pos = [i for i in range(len(tip) - len(kmer) + 1) if tip.startswith(kmer, i)]
    spacings = Counter(b - a for a, b in zip(pos, pos[1:], strict=False))
    return spacings.most_common(1)[0][0] if spacings else None


def _oriented_sequence(graph: BidirectedGraph, sym: str) -> str:
    seq = graph.segments[sym[:-1]].sequence or ""
    return seq if sym.endswith("+") else reverse_complement(seq)


def build_end_catalog(
    graph: BidirectedGraph,
    ledger_map: Mapping[str, Any],
    agp_records: Iterable[Any],
    orphan_caps: Mapping[str, tuple[str, str]],
    core_coverage: float,
    paf_hits: Sequence[PafHit] | None = None,
    policy: Policy = DEFAULT_POLICY,
    end_evidence: PhasedEndEvidence | None = None,
) -> tuple[list[EndCatalogEntry], list[AmbiguityGroup]]:
    """Assess every end of the nominated chromosome objects.

    Args:
        graph: Primary unitig graph with sequences loaded.
        ledger_map: Port symbol to `EndEvidence`.
        agp_records: Emitted AGP records after SUPER nomination.
        orphan_caps: Segment name to (telomere port, inward port), from `identify_orphan_caps`.
        core_coverage: Median core depth used to call core-depth unitigs.
        paf_hits: Optional long-read alignments to the unitigs.
        policy: Thresholds shared with `finish`.

    Returns:
        Catalog entries ordered by object and side, and the ambiguity groups.
    """
    min_core = policy.stoichiometry.min_stem_core_ratio.value * core_coverage

    def is_core(name: str) -> bool:
        seg = graph.segments.get(name)
        return seg is not None and core_coverage > 0 and (seg.depth or 0.0) >= min_core

    def has_telomere(sym: str) -> bool:
        ev = ledger_map.get(sym)
        return bool(ev is not None and ev.has_terminal_telomere)

    def has_rdna(name: str) -> bool:
        return any(getattr(ledger_map.get(f"{name}{s}"), "rdna_hit_count", 0) > 0 for s in "+-")

    by_object: dict[str, list[Any]] = {}
    for rec in agp_records:
        if rec.object_id.startswith(CHROMOSOME_PREFIX):
            by_object.setdefault(rec.object_id, []).append(rec)

    entries: list[EndCatalogEntry] = []
    for object_id in sorted(by_object, key=lambda o: int(o[len(CHROMOSOME_PREFIX) :])):
        for side in ("5p", "3p"):
            items = _outward_items(by_object[object_id], side, graph)
            if items is not None:
                entries.append(_walk_end(object_id, side, items, graph, is_core, has_telomere))

    certain_units = {p[:-1] for e in entries for p in e.certain_ports}
    pool = {name: tel for name, (tel, _) in orphan_caps.items() if name not in certain_units}
    window = policy.ont.min_anchor_bp.value * 2
    hits_by_target: dict[str, list[PafHit]] = {}
    hits_by_read: dict[str, list[PafHit]] = {}
    for hit in paf_hits or ():
        hits_by_target.setdefault(hit.tname, []).append(hit)
        hits_by_read.setdefault(hit.qname, []).append(hit)
    engine = OntBridgingEngine(graph, ledger_map, policy)

    families_by_id: dict[str, EndEvidenceFamily] = {}
    placements_by_record: dict[str, list[ConsensusPlacement]] = {}
    records_by_cluster: dict[str, set[str]] = {}
    clusters_by_molecule: dict[str, set[str]] = {}
    if end_evidence is not None:
        families_by_id = {f.cluster_id: f for f in end_evidence.families}
        for fam in end_evidence.families:
            for member in fam.members:
                clusters_by_molecule.setdefault(member.molecule_id, set()).add(fam.cluster_id)
        for pl in end_evidence.consensus_placements:
            if pl.source != "primary":
                continue
            placements_by_record.setdefault(pl.record, []).append(pl)
            records_by_cluster.setdefault(pl.cluster_id, set()).add(pl.record)

    for entry in entries:
        if entry.horizon_reason is None:
            continue
        evidence: dict[str, list[str]] = {}
        candidate_families: dict[str, set[str]] = {}
        candidate_mols: dict[str, set[str]] = {}
        block: set[str] = set(entry.block_unitigs)

        # Graph: everything reachable outward without re-entering a certain path or the core.
        seen: set[str] = set()
        queue = deque([entry.horizon_port])
        while queue:
            for nb in graph.neighbors(_port(queue.popleft())):
                sym = nb.to_symbol()
                if nb.name in certain_units or sym in seen:
                    continue
                seen.add(sym)
                block.add(nb.name)
                queue.append(sym)
                if has_telomere(sym):
                    evidence.setdefault(sym, []).append("graph")

        anchored: list[tuple[PafHit, int, int]] = []
        if hits_by_target:
            for name, tel in sorted(pool.items()):
                subset = hits_by_target.get(entry.horizon_port[:-1], []) + hits_by_target.get(name, [])
                spanning = engine.collect_spanning_reads(subset, entry.horizon_port, tel)
                if spanning:
                    evidence.setdefault(tel, []).append(f"reads:{len(spanning)}")
                    candidate_mols.setdefault(tel, set()).update(sr.qname for sr in spanning)
            anchored = _horizon_reads(entry, hits_by_target, policy)
            anchor_reads: dict[str, set[str]] = {}
            for hit, outward_from, outward_sign in anchored:
                for other in hits_by_read[hit.qname]:
                    if other.tname not in pool or other.mapq < policy.ont.min_catalog_link_mapq.value:
                        continue
                    if (other.qstart + other.qend) / 2 * outward_sign > outward_from * outward_sign:
                        anchor_reads.setdefault(pool[other.tname], set()).add(hit.qname)
            for tel, reads in sorted(anchor_reads.items()):
                evidence.setdefault(tel, []).append(f"read-anchor:{len(reads)}")
                candidate_mols.setdefault(tel, set()).update(reads)

        oriented = _oriented_sequence(graph, entry.horizon_port)
        tip = tip_satellite_kmer(oriented, window, policy)
        if tip is not None:
            kmer, rc_kmer = tip[0], reverse_complement(tip[0])
            entry.tip_repeat_kmer = kmer
            entry.tip_repeat_period_bp = _tip_period(oriented, kmer, window)
            floor = policy.ont.min_spanning_reads.value
            for name, tel in sorted(pool.items()):
                seq = graph.segments[name].sequence or ""
                if seq.count(kmer) + seq.count(rc_kmer) >= floor:
                    evidence.setdefault(tel, []).append("satellite-kmer")
            for name, seg in graph.segments.items():
                if name in certain_units or name in orphan_caps or not seg.sequence:
                    continue
                if graph.degree(SegmentPort(name, Strand.FORWARD)) or graph.degree(SegmentPort(name, Strand.REVERSE)):
                    continue
                if seg.sequence.count(kmer) + seg.sequence.count(rc_kmer) >= SATELLITE_MIN_TIP_REPEATS:
                    block.add(name)

        if entry.horizon_reason == HorizonReason.RDNA:
            for name, tel in sorted(pool.items()):
                if has_rdna(name):
                    evidence.setdefault(tel, []).append("rdna")

        if end_evidence is not None:
            entry_clusters: set[str] = set()
            for port_sym in entry.certain_ports:
                for pl in placements_by_record.get(port_sym[:-1], ()):
                    entry_clusters.add(pl.cluster_id)
            for hit, _, _ in anchored:
                entry_clusters.update(clusters_by_molecule.get(hit.qname, ()))

            for cid in sorted(entry_clusters):
                for name, tel in sorted(pool.items()):
                    if name in records_by_cluster.get(cid, ()):
                        evidence.setdefault(tel, [])

            for tel in list(evidence):
                for pl in placements_by_record.get(tel[:-1], ()):
                    cid = pl.cluster_id
                    ev_tag = f"family:{cid}"
                    if ev_tag not in evidence[tel]:
                        evidence[tel].append(ev_tag)
                    candidate_families.setdefault(tel, set()).add(cid)
                    fam = families_by_id.get(cid)
                    if fam is not None:
                        candidate_mols.setdefault(tel, set()).update(m.molecule_id for m in fam.members)

        entry.candidates = [
            CandidateEnding(
                p,
                tuple(ev),
                family_ids=tuple(sorted(candidate_families.get(p, ()))),
                molecule_ids=tuple(sorted(candidate_mols.get(p, ()))),
            )
            for p, ev in sorted(evidence.items())
        ]
        entry.block_unitigs = tuple(sorted(block - {c.port[:-1] for c in entry.candidates}))

    groups = _group(entries, pool)
    return entries, groups


def _walk_end(
    object_id: str,
    side: str,
    items: list[_Item],
    graph: BidirectedGraph,
    is_core: Any,
    has_telomere: Any,
) -> EndCatalogEntry:
    """Split the outward component list at the first join the data do not establish."""
    certain = [items[0].port]
    joins: list[str] = []
    reason: HorizonReason | None = None
    stop = len(items)
    prev_port = items[0].port
    gap_join: str | None = None
    for i in range(1, len(items)):
        cur = items[i]
        if cur.kind != "W":
            gap_join = _CERTAIN_GAP_JOIN.get(cur.linkage or "")
            if gap_join is None:
                reason = _LINKAGE_REASON.get(cur.linkage or "", HorizonReason.GAP)
                stop = i
                break
            continue
        assert prev_port is not None and cur.port is not None
        if gap_join is not None:
            joins.append(gap_join)
            gap_join = None
        elif prev_port[:-1] not in graph.segments or cur.port[:-1] not in graph.segments:
            joins.append("read")
        elif cur.port in {n.to_symbol() for n in graph.neighbors(_port(prev_port))}:
            joins.append("graph")
        else:
            reason = HorizonReason.NO_GRAPH_EDGE
            stop = i
            break
        certain.append(cur.port)
        prev_port = cur.port

    unitig_ports = [p for p in certain if p[:-1] in graph.segments]
    horizon = unitig_ports[-1]
    beyond = [it.port for it in items[stop:] if it.port]
    finish_choice = beyond[-1] if beyond and has_telomere(beyond[-1]) else None
    if reason is None:
        if not has_telomere(horizon):
            reason = HorizonReason.NO_TELOMERE
        elif is_core(horizon[:-1]):
            reason = HorizonReason.CORE_TELOMERE
    return EndCatalogEntry(
        object_id=object_id,
        side=side,
        certain_ports=tuple(certain),
        certain_joins=tuple(joins),
        horizon_port=horizon,
        horizon_reason=reason,
        finish_choice=finish_choice,
        finish_linkage=items[stop].linkage if stop < len(items) and items[stop].kind != "W" else None,
        block_unitigs=tuple(p[:-1] for p in beyond if p != finish_choice and p[:-1] in graph.segments),
    )


def _horizon_reads(
    entry: EndCatalogEntry, hits_by_target: Mapping[str, list[PafHit]], policy: Policy
) -> list[tuple[PafHit, int, int]]:
    """Reads anchored at the horizon boundary, with where and which way they leave it.

    Returns ``(hit, read coordinate at the boundary, +1 or -1 outward along the read)``
    and records the anchored read count and longest extension past the boundary.
    """
    port = _port(entry.horizon_port)
    anchored: list[tuple[PafHit, int, int]] = []
    overhangs: list[int] = []
    for hit in hits_by_target.get(port.name, []):
        if hit.alen < policy.ont.min_anchor_bp.value or hit.mapq < policy.ont.min_mapq.value:
            continue
        slack = _port_boundary_slack(hit, True, port.strand)
        if slack > policy.ont.max_boundary_slack_bp.value:
            continue
        same = (hit.strand == "+") == (port.strand == Strand.FORWARD)
        overhangs.append(max(0, (hit.qlen - hit.qend if same else hit.qstart) - slack))
        anchored.append((hit, hit.qend, 1) if same else (hit, hit.qstart, -1))
    entry.reads_at_horizon = len({hit.qname for hit, _, _ in anchored})
    entry.max_read_overhang_bp = max(overhangs) if overhangs else None
    return anchored


def _group(entries: list[EndCatalogEntry], pool: Mapping[str, str]) -> list[AmbiguityGroup]:
    """Union ends that share a candidate cap or an unresolved unitig; list unclaimed caps."""
    parent: dict[str, str] = {}

    def find(x: str) -> str:
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    open_entries = [e for e in entries if e.horizon_reason is not None]
    for e in open_entries:
        find(e.end_id)
        for name in [c.port[:-1] for c in e.candidates] + list(e.block_unitigs):
            parent[find(f"utg:{name}")] = find(e.end_id)

    members: dict[str, list[EndCatalogEntry]] = {}
    for e in open_entries:
        members.setdefault(find(e.end_id), []).append(e)

    groups: list[AmbiguityGroup] = []
    for idx, group_entries in enumerate(members.values(), start=1):
        gid = f"G{idx}"
        caps = sorted({c.port for e in group_entries for c in e.candidates})
        unitig_ends: dict[str, int] = Counter(
            name for e in group_entries for name in {c.port[:-1] for c in e.candidates} | set(e.block_unitigs)
        )
        for e in group_entries:
            e.group_id = gid
            if len(group_entries) > 1:
                e.category = EndCategory.SHARED
            elif e.candidates:
                e.category = EndCategory.BOUNDED
            else:
                e.category = EndCategory.OPEN
        groups.append(
            AmbiguityGroup(
                group_id=gid,
                ends=tuple(e.end_id for e in group_entries),
                candidate_endings=tuple(caps),
                shared_unitigs=tuple(sorted(n for n, c in unitig_ends.items() if c > 1)),
            )
        )

    claimed = {c.port[:-1] for e in open_entries for c in e.candidates}
    for name in sorted(set(pool) - claimed):
        groups.append(AmbiguityGroup(f"G{len(groups) + 1}", (), (pool[name],), ()))
    return groups


END_CATALOG_COLUMNS = (
    "end_id",
    "object_id",
    "side",
    "category",
    "horizon_port",
    "horizon_reason",
    "certain_ports",
    "certain_joins",
    "candidate_endings",
    "n_candidates",
    "group_id",
    "finish_choice",
    "finish_linkage",
    "block_unitigs",
    "tip_repeat_kmer",
    "tip_repeat_period_bp",
    "reads_at_horizon",
    "max_read_overhang_bp",
    "candidate_families",
    "candidate_molecules",
)

GROUP_COLUMNS = ("group_id", "n_ends", "ends", "candidate_endings", "shared_unitigs", "fusion_interpretations")


def _cell(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, tuple | list):
        return ",".join(str(v) for v in value)
    return str(value)


def emit_end_catalog(entries: Sequence[EndCatalogEntry], groups: Sequence[AmbiguityGroup], out_dir: Path) -> None:
    """Write `end_catalog.tsv` and `ambiguity_groups.tsv` into ``out_dir``."""
    with (out_dir / "end_catalog.tsv").open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(END_CATALOG_COLUMNS)
        for e in entries:
            writer.writerow(
                [
                    e.end_id,
                    e.object_id,
                    e.side,
                    e.category.value,
                    e.horizon_port,
                    _cell(e.horizon_reason and e.horizon_reason.value),
                    _cell(e.certain_ports),
                    _cell(e.certain_joins),
                    ";".join(c.format() for c in e.candidates),
                    len(e.candidates),
                    _cell(e.group_id),
                    _cell(e.finish_choice),
                    _cell(e.finish_linkage),
                    _cell(e.block_unitigs),
                    _cell(e.tip_repeat_kmer),
                    _cell(e.tip_repeat_period_bp),
                    _cell(e.reads_at_horizon),
                    _cell(e.max_read_overhang_bp),
                    ";".join(",".join(c.family_ids) for c in e.candidates)
                    if any(c.family_ids for c in e.candidates)
                    else "",
                    ";".join(",".join(c.molecule_ids) for c in e.candidates)
                    if any(c.molecule_ids for c in e.candidates)
                    else "",
                ]
            )
    with (out_dir / "ambiguity_groups.tsv").open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(GROUP_COLUMNS)
        for g in groups:
            writer.writerow(
                [
                    g.group_id,
                    len(g.ends),
                    _cell(g.ends),
                    _cell(g.candidate_endings),
                    _cell(g.shared_unitigs),
                    _cell(g.fusion_interpretations),
                ]
            )


__all__ = [
    "AmbiguityGroup",
    "CandidateEnding",
    "ConsensusPlacement",
    "EndCatalogEntry",
    "EndCategory",
    "EndEvidenceFamily",
    "EndEvidenceMember",
    "HaplotypePattern",
    "HorizonReason",
    "PhasedEndEvidence",
    "build_end_catalog",
    "emit_end_catalog",
    "load_end_evidence",
]
