"""Core domain data models for bidirected assembly graphs.

Provides foundational representations for strands, segment ports, segments,
and dovetail overlap edges with Watson-Crick duality invariance.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Any


class _StrandCanonicalizeDescriptor:
    """Descriptor supporting canonicalize on both Strand class and instances."""

    def __get__(self, instance: Any, owner: type[Strand]) -> Any:
        def _canonicalize(val: str | Strand | None = None) -> Strand:
            if val is not None:
                return owner.from_symbol(val)
            if instance is not None:
                return instance
            return owner.FORWARD

        return _canonicalize


class Strand(str, Enum):
    """Sequence or port orientation strand (+ or -)."""

    FORWARD = "+"
    REVERSE = "-"
    PLUS = "+"
    MINUS = "-"
    POS = "+"
    NEG = "-"

    canonicalize = _StrandCanonicalizeDescriptor()

    def invert(self) -> Strand:
        """Return the inverted (opposite) strand."""
        return Strand.REVERSE if self == Strand.FORWARD else Strand.FORWARD

    def __invert__(self) -> Strand:
        """Bitwise not (~strand) returns the inverted strand."""
        return self.invert()

    def to_symbol(self) -> str:
        """Return the single-character strand symbol ('+' or '-')."""
        return self.value

    @classmethod
    def from_symbol(cls, symbol: str | Strand) -> Strand:
        """Resolve a strand symbol or instance to a canonical Strand."""
        if isinstance(symbol, cls):
            return symbol
        if symbol in ("+", "1", "fwd", "forward", "plus", "POS"):
            return cls.FORWARD
        if symbol in ("-", "-1", "rev", "reverse", "minus", "NEG"):
            return cls.REVERSE
        raise ValueError(f"Invalid strand symbol: {symbol!r}")


@dataclass(slots=True, frozen=True)
class SegmentPort:
    """An oriented port of a segment in a bidirected graph."""

    name: str
    strand: Strand

    def __post_init__(self) -> None:
        if isinstance(self.strand, str) and not isinstance(self.strand, Strand):
            object.__setattr__(self, "strand", Strand.from_symbol(self.strand))

    def invert(self) -> SegmentPort:
        """Return the reverse-complement port."""
        return SegmentPort(name=self.name, strand=self.strand.invert())

    def __invert__(self) -> SegmentPort:
        return self.invert()

    @property
    def segment_id(self) -> str:
        """Alias for segment name."""
        return self.name

    def to_symbol(self) -> str:
        """Return port representation in name+ or name- notation."""
        return f"{self.name}{self.strand.to_symbol()}"


@dataclass(slots=True)
class Segment:
    """A unitig or sequence segment in the assembly graph."""

    name: str
    length: int
    depth: float = 0.0
    sequence: str | None = None

    def __len__(self) -> int:
        return self.length


@dataclass(slots=True, frozen=True)
class DovetailEdge:
    """A directed dovetail overlap edge between two segment ports.

    Follows GFA1 link semantics where an edge from u(strand) to v(strand)
    describes an overlap between the 3' suffix of u and the 5' prefix of v.
    """

    from_port: SegmentPort
    to_port: SegmentPort
    cigar: str
    overlap_len: int
    l1: int = 0
    l2: int = 0

    @property
    def dual(self) -> DovetailEdge:
        """Return the reverse-complement dual edge in the bidirected graph.

        Dual of (u -> v, l1, l2) is (~v -> ~u, l2, l1).
        Satisfies edge.dual.dual == edge (involution invariant).
        """
        return DovetailEdge(
            from_port=self.to_port.invert(),
            to_port=self.from_port.invert(),
            cigar=self.cigar,
            overlap_len=self.overlap_len,
            l1=self.l2,
            l2=self.l1,
        )

    def __invert__(self) -> DovetailEdge:
        return self.dual

    def __call__(self) -> DovetailEdge:
        """Permit dual to be accessed either as property (edge.dual) or method (edge.dual())."""
        return self
