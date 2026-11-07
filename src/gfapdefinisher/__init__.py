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
    "DEFAULT_POLICY",
    "EndState",
    "EndStateRecord",
    "GapDistribution",
    "IntervalQuirk",
    "JunctionClass",
    "JunctionRecord",
    "OntBridgingEngine",
    "PafHit",
    "PdeJunctionCall",
    "Policy",
    "PortClassification",
    "Provenance",
    "QuirkRecord",
    "RdnaCassetteResolver",
    "RdnaResolution",
    "SatelliteGapResolution",
    "SpanningRead",
    "TerminalFeature",
    "Topology",
    "classify_pde_junction",
    "classify_port",
    "hello",
    "parse_paf_hits",
]
