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


def hello() -> str:
    return "Hello from gfapdefinisher!"


__all__ = [
    "PdeJunctionCall",
    "PortClassification",
    "Provenance",
    "TerminalFeature",
    "Topology",
    "classify_pde_junction",
    "classify_port",
    "hello",
]
