"""gfapdefinisher: T2T assembly finishing and Programmed DNA Elimination harness."""

from gfapdefinisher.classifier import (
    PdeJunctionCall,
    PortClassification,
    Provenance,
    TerminalFeature,
    Topology,
    classify_pde_junction,
    classify_port,
)
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
    emit_agp,
    emit_fasta,
    emit_scaffold,
    parse_agp,
)
from gfapdefinisher.policy import (
    DEFAULT_POLICY,
    Policy,
)
from gfapdefinisher.quirks import (
    EndState,
    EndStateRecord,
    IntervalQuirk,
    JunctionClass,
    JunctionRecord,
    QuirkRecord,
)
from gfapdefinisher.rdna import (
    RdnaCassetteResolver,
    RdnaResolution,
)
from gfapdefinisher.resolver import (
    GapDistribution,
    OntBridgingEngine,
    PafHit,
    SatelliteGapResolution,
    SpanningRead,
    parse_paf_hits,
)


def hello() -> str:
    return "Hello from gfapdefinisher!"


__all__ = [
    "AgpGapRow",
    "AgpRecord",
    "DEFAULT_POLICY",
    "EndState",
    "EndStateRecord",
    "FastaAgpEmitter",
    "GapDistribution",
    "IntervalQuirk",
    "JunctionClass",
    "JunctionRecord",
    "OntBridgingEngine",
    "PathComponent",
    "PafHit",
    "PdeJunctionCall",
    "Policy",
    "PortClassification",
    "Provenance",
    "QuirkRecord",
    "RdnaCassetteResolver",
    "RdnaResolution",
    "SatelliteGapResolution",
    "ScaffoldPlan",
    "ScaffoldStep",
    "SequenceProvenance",
    "SpanningRead",
    "TerminalFeature",
    "Topology",
    "assemble_path_sequence",
    "assemble_scaffold_sequence",
    "classify_pde_junction",
    "classify_port",
    "emit_agp",
    "emit_fasta",
    "emit_scaffold",
    "hello",
    "parse_agp",
    "parse_paf_hits",
]
