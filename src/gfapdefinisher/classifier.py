"""Port classification and Programmed DNA Elimination (PDE) junction motif detector.

Implements 3-axis port classification (Provenance, Topology, TerminalFeature)
and PDE junction motif detection for bidirected assembly graphs.
Adheres strictly to the single source of truth in docs/biology/CLASSIFIERS.md
and biological specifications BIO-001 through BIO-006:
- BIO-001: Authentic telomeric cap reachability and completion.
- BIO-002: Stoichiometric depth contrast over arbitrary length cutoffs.
- BIO-003: Graph-inferred diminution boundaries without external coordinate dependencies.
- BIO-004: Repeat representation preservation for variable satellites and complex GRS.
- BIO-005: Lossless overlap and provenance invariance.
- BIO-006: Generalization across nematode genomes without benchmark over-tuning.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Any

from gfapdefinisher.graph import BidirectedGraph
from gfapdefinisher.ledger import EndEvidence
from gfapdefinisher.models import SegmentPort, Strand
from gfapdefinisher.policy import DEFAULT_POLICY, Policy


class Provenance(str, Enum):
    """Taxonomic and genomic provenance origin of a segment."""

    NUCLEAR = "nuclear"
    ORGANELLE = "organelle"
    CONTAMINANT = "contaminant"
    UNKNOWN = "unknown"


class Topology(str, Enum):
    """Outgoing graph topology degree from an oriented segment port."""

    TIP = "tip"
    COLLINEAR = "collinear"
    FORK = "fork"
    COMPLEX_BRANCH = "complex_branch"


class TerminalFeature(str, Enum):
    """Biological sequence composition feature observed at a segment terminus."""

    TELOMERIC_CAP = "telomeric_cap"
    RDNA_CLUSTER = "rdna_cluster"
    NONE = "none"


@dataclass(slots=True, frozen=True)
class PortClassification:
    """Three-axis orthogonal classification for an oriented segment port."""

    port: SegmentPort
    provenance: Provenance
    topology: Topology
    terminal_feature: TerminalFeature


@dataclass(slots=True, frozen=True)
class PdeJunctionCall:
    """Decision call evaluating whether an oriented port is a PDE junction.

    See docs/biology/CLASSIFIERS.md Section 2 for the formal specification.
    """

    port: SegmentPort
    is_pde_junction: bool
    somatic_port: SegmentPort | None = None
    germline_port: SegmentPort | None = None
    depth_ratio: float = 0.0


def _resolve_port(port: SegmentPort | str) -> SegmentPort:
    """Resolve a port representation to a canonical SegmentPort."""
    if isinstance(port, SegmentPort):
        return port
    if port.endswith("+"):
        return SegmentPort(name=port[:-1], strand=Strand.FORWARD)
    if port.endswith("-"):
        return SegmentPort(name=port[:-1], strand=Strand.REVERSE)
    return SegmentPort(name=port, strand=Strand.FORWARD)


def _lookup_evidence(
    ledger_by_port: dict[str, EndEvidence] | dict[Any, EndEvidence],
    port: SegmentPort,
) -> EndEvidence | None:
    """Look up EndEvidence for a port across supported indexing keys."""
    if port in ledger_by_port:
        return ledger_by_port[port]

    sym = port.to_symbol()
    if sym in ledger_by_port:
        return ledger_by_port[sym]

    end_str = "5p" if port.strand == Strand.MINUS else "3p"
    colon_key = f"{port.name}:{end_str}"
    if colon_key in ledger_by_port:
        return ledger_by_port[colon_key]

    strand_key = f"{port.name}:{port.strand.to_symbol()}"
    if strand_key in ledger_by_port:
        return ledger_by_port[strand_key]

    return None


def _is_telomere_capped(
    port: SegmentPort,
    graph: BidirectedGraph,
    ledger_by_port: dict[str, EndEvidence] | dict[Any, EndEvidence],
) -> bool:
    """Check whether a port terminates in an authentic telomeric cap tip (BIO-001)."""
    ev = _lookup_evidence(ledger_by_port, port)
    is_dead_end = graph.degree(port) == 0 if port.name in graph.segments else (ev.is_tip if ev is not None else False)
    has_cap = ev.has_terminal_telomere if ev is not None else False
    return is_dead_end and has_cap


def _reaches_somatic_telomere_cap(
    start_port: SegmentPort,
    graph: BidirectedGraph,
    ledger_by_port: dict[str, EndEvidence] | dict[Any, EndEvidence],
    stem_depth: float,
    max_germline_ratio: float,
    max_hops: int,
) -> bool:
    """Evaluate if start_port reaches a somatic telomeric cap within max_hops traversal.

    Governed by policy.pde.max_somatic_cap_hops and docs/biology/CLASSIFIERS.md Section 2.
    """
    if _is_telomere_capped(start_port, graph, ledger_by_port):
        return True

    if max_hops <= 0:
        return False

    queue: list[tuple[SegmentPort, int]] = [(start_port, 0)]
    visited: set[str] = {start_port.to_symbol()}

    while queue:
        curr_port, hops = queue.pop(0)
        if hops >= max_hops:
            continue

        for edge in graph.get_outgoing(curr_port):
            nxt_port = edge.to_port
            nxt_sym = nxt_port.to_symbol()
            if nxt_sym in visited:
                continue
            visited.add(nxt_sym)

            nxt_seg = graph.segments.get(nxt_port.name)
            nxt_ev = _lookup_evidence(ledger_by_port, nxt_port)
            if nxt_seg is not None and nxt_seg.depth > 0:
                nxt_depth = nxt_seg.depth
            elif nxt_ev is not None and nxt_ev.depth > 0:
                nxt_depth = nxt_ev.depth
            else:
                nxt_depth = 0.0

            if (nxt_depth / stem_depth) > max_germline_ratio:
                if _is_telomere_capped(nxt_port, graph, ledger_by_port):
                    return True
                queue.append((nxt_port, hops + 1))

    return False


def classify_port(
    port: SegmentPort | str,
    graph: BidirectedGraph,
    evidence: EndEvidence,
) -> PortClassification:
    """Assign orthogonal 3-axis labels (Provenance, Topology, TerminalFeature) to a port.

    Canonical taxonomy defined in docs/biology/CLASSIFIERS.md Section 1.

    Args:
        port: Oriented segment port being classified.
        graph: Bidirected assembly graph context.
        evidence: Tier B end-evidence ledger entry for the port.

    Returns:
        PortClassification containing the 3-axis taxonomic and topological call.
    """
    resolved = _resolve_port(port)

    # Axis 1: Provenance
    tiara_cls = (evidence.tiara_class or "").lower().strip()
    if evidence.usco_hit_count > 0 or tiara_cls == "eukarya":
        provenance = Provenance.NUCLEAR
    elif tiara_cls in ("organelle", "mitochondrion", "plastid"):
        provenance = Provenance.ORGANELLE
    elif tiara_cls in ("bacteria", "archaea", "contaminant"):
        provenance = Provenance.CONTAMINANT
    else:
        provenance = Provenance.UNKNOWN

    # Axis 2: Topology
    deg = graph.degree(resolved) if resolved.name in graph.segments else evidence.degree
    if deg == 0:
        topology = Topology.TIP
    elif deg == 1:
        topology = Topology.COLLINEAR
    elif deg == 2:
        topology = Topology.FORK
    else:
        topology = Topology.COMPLEX_BRANCH

    # Axis 3: TerminalFeature (BIO-001)
    if evidence.has_terminal_telomere:
        terminal_feature = TerminalFeature.TELOMERIC_CAP
    elif evidence.has_terminal_rdna:
        terminal_feature = TerminalFeature.RDNA_CLUSTER
    else:
        terminal_feature = TerminalFeature.NONE

    return PortClassification(
        port=resolved,
        provenance=provenance,
        topology=topology,
        terminal_feature=terminal_feature,
    )


def classify_pde_junction(
    port: SegmentPort | str,
    graph: BidirectedGraph,
    ledger_by_port: dict[str, EndEvidence] | dict[Any, EndEvidence],
    policy: Policy = DEFAULT_POLICY,
) -> PdeJunctionCall:
    """Evaluate whether an oriented segment port is a Programmed DNA Elimination (PDE) junction.

    See docs/biology/CLASSIFIERS.md Section 2 for the formal biological specification.
    A PDE junction is characterized by the reachability and stoichiometry motif (BIO-002, BIO-003):
    1. A somatic stem with robust read depth.
    2. Bifurcation into a somatic branch reaching an authentic telomeric cap within <= max_somatic_cap_hops
       with read depth near/above somatic levels (depth ratio > max_germline_depth_ratio).
    3. Concurrently bifurcating into a germline-retained sequence (GRS) branch with stoichiometric
       depth ratio in [min_germline_depth_ratio, max_germline_depth_ratio] (BIO-002, BIO-004).

    Args:
        port: Stem port to evaluate.
        graph: Bidirected assembly graph.
        ledger_by_port: Pre-computed end-evidence mapping.
        policy: Operational parameter policy governing stoichiometry and reachability.

    Returns:
        PdeJunctionCall with somatic_port and germline_port references when identified.
    """
    resolved_stem = _resolve_port(port)

    # Stem must bifurcate or form a complex branch (deg >= 2)
    stem_deg = graph.degree(resolved_stem) if resolved_stem.name in graph.segments else 0
    if stem_deg < 2:
        return PdeJunctionCall(port=resolved_stem, is_pde_junction=False)

    stem_seg = graph.segments.get(resolved_stem.name)
    stem_ev = _lookup_evidence(ledger_by_port, resolved_stem)
    if stem_seg is not None and stem_seg.depth > 0:
        stem_depth = stem_seg.depth
    elif stem_ev is not None and stem_ev.depth > 0:
        stem_depth = stem_ev.depth
    else:
        stem_depth = 0.0

    if stem_depth <= 0.0:
        return PdeJunctionCall(port=resolved_stem, is_pde_junction=False)

    min_germline_ratio = policy.stoichiometry.min_germline_depth_ratio.value
    max_germline_ratio = policy.stoichiometry.max_germline_depth_ratio.value
    max_hops = policy.pde.max_somatic_cap_hops.value

    outgoing_edges = graph.get_outgoing(resolved_stem)

    somatic_candidates: list[tuple[SegmentPort, float]] = []
    germline_candidates: list[tuple[SegmentPort, float]] = []

    for edge in outgoing_edges:
        branch_port = edge.to_port
        branch_seg = graph.segments.get(branch_port.name)
        branch_ev = _lookup_evidence(ledger_by_port, branch_port)

        if branch_seg is not None and branch_seg.depth > 0:
            branch_depth = branch_seg.depth
        elif branch_ev is not None and branch_ev.depth > 0:
            branch_depth = branch_ev.depth
        else:
            branch_depth = 0.0

        ratio = branch_depth / stem_depth

        # Somatic candidate: high depth ratio and reaches telomeric cap within max_hops
        if ratio > max_germline_ratio:
            if _reaches_somatic_telomere_cap(
                branch_port,
                graph,
                ledger_by_port,
                stem_depth,
                max_germline_ratio,
                max_hops,
            ):
                somatic_candidates.append((branch_port, ratio))

        # Germline candidate: stoichiometric depth ratio in [min_ratio, max_ratio] (BIO-002, BIO-004)
        if min_germline_ratio <= ratio <= max_germline_ratio:
            germline_candidates.append((branch_port, ratio))

    if somatic_candidates and germline_candidates:
        # Deterministic candidate selection
        somatic_candidates.sort(key=lambda item: (-item[1], item[0].name, item[0].strand.value))
        germline_candidates.sort(key=lambda item: (item[0].name, item[0].strand.value))

        best_somatic_port, _ = somatic_candidates[0]
        best_germline_port, best_germline_ratio = germline_candidates[0]

        return PdeJunctionCall(
            port=resolved_stem,
            is_pde_junction=True,
            somatic_port=best_somatic_port,
            germline_port=best_germline_port,
            depth_ratio=best_germline_ratio,
        )

    return PdeJunctionCall(port=resolved_stem, is_pde_junction=False)
