"""Unit tests for core data models (Strand, SegmentPort, Segment, DovetailEdge)."""

from dataclasses import FrozenInstanceError

import pytest

from gfapdefinisher.models import DovetailEdge, Segment, SegmentPort, Strand


def test_strand_values_and_inversion():
    """Verify Strand enum values, string representation, and inversion involution."""
    assert Strand.FORWARD == "+"
    assert Strand.REVERSE == "-"
    assert Strand.FORWARD.to_symbol() == "+"
    assert Strand.REVERSE.to_symbol() == "-"

    # Invert method
    assert Strand.FORWARD.invert() == Strand.REVERSE
    assert Strand.REVERSE.invert() == Strand.FORWARD

    # Invert operator
    assert ~Strand.FORWARD == Strand.REVERSE
    assert ~Strand.REVERSE == Strand.FORWARD

    # Double invert invariant
    assert Strand.FORWARD.invert().invert() == Strand.FORWARD
    assert ~(~Strand.REVERSE) == Strand.REVERSE


def test_strand_aliases_and_from_symbol():
    """Verify alias resolution and symbol parsing."""
    assert Strand.PLUS == Strand.FORWARD
    assert Strand.MINUS == Strand.REVERSE
    assert Strand.POS == Strand.FORWARD
    assert Strand.NEG == Strand.REVERSE

    assert Strand.from_symbol("+") == Strand.FORWARD
    assert Strand.from_symbol("-") == Strand.REVERSE
    assert Strand.from_symbol("fwd") == Strand.FORWARD
    assert Strand.from_symbol("rev") == Strand.REVERSE
    assert Strand.from_symbol(Strand.FORWARD) == Strand.FORWARD

    with pytest.raises(ValueError, match="Invalid strand symbol"):
        Strand.from_symbol("?")


def test_strand_canonicalize():
    """Verify canonicalize works as classmethod and instance method."""
    assert Strand.canonicalize("+") == Strand.FORWARD
    assert Strand.canonicalize("-") == Strand.REVERSE
    assert Strand.canonicalize() == Strand.FORWARD
    assert Strand.FORWARD.canonicalize() == Strand.FORWARD
    assert Strand.REVERSE.canonicalize() == Strand.REVERSE


def test_segment_port_creation_and_inversion():
    """Verify SegmentPort attributes, conversions, and inversion."""
    p_fwd = SegmentPort(name="utg001", strand=Strand.FORWARD)
    assert p_fwd.name == "utg001"
    assert p_fwd.strand == Strand.FORWARD
    assert p_fwd.to_symbol() == "utg001+"

    # String strand auto-conversion
    p_str = SegmentPort(name="utg002", strand="-")
    assert p_str.strand == Strand.REVERSE
    assert p_str.to_symbol() == "utg002-"

    # Inversion
    inv = p_fwd.invert()
    assert inv.name == "utg001"
    assert inv.strand == Strand.REVERSE
    assert ~p_fwd == inv

    # Double invert invariant
    assert p_fwd.invert().invert() == p_fwd
    assert ~(~p_fwd) == p_fwd


def test_segment_port_immutability_and_hashing():
    """Verify SegmentPort is frozen and hashable."""
    p1 = SegmentPort(name="utg001", strand=Strand.FORWARD)
    p2 = SegmentPort(name="utg001", strand=Strand.FORWARD)
    p3 = SegmentPort(name="utg001", strand=Strand.REVERSE)

    assert p1 == p2
    assert hash(p1) == hash(p2)
    assert p1 != p3

    # Set storage
    port_set = {p1, p2, p3}
    assert len(port_set) == 2

    with pytest.raises(FrozenInstanceError):
        p1.name = "changed"  # type: ignore[misc]


def test_segment_model():
    """Verify Segment dataclass attributes and length reporting."""
    seg = Segment(name="utg001", length=15000, depth=32.5, sequence="ACGT" * 10)
    assert seg.name == "utg001"
    assert seg.length == 15000
    assert seg.depth == 32.5
    assert len(seg) == 15000
    assert seg.sequence is not None and len(seg.sequence) == 40

    # Defaults
    seg_default = Segment(name="utg002", length=500)
    assert seg_default.depth == 0.0
    assert seg_default.sequence is None


def test_dovetail_edge_and_dual_invariant():
    """Verify DovetailEdge creation, dual calculation, and involution invariant."""
    p1 = SegmentPort(name="utg001", strand=Strand.FORWARD)
    p2 = SegmentPort(name="utg002", strand=Strand.FORWARD)

    edge = DovetailEdge(
        from_port=p1,
        to_port=p2,
        cigar="15643M",
        overlap_len=15643,
        l1=20000,
        l2=30000,
    )

    # Edge attributes
    assert edge.from_port == p1
    assert edge.to_port == p2
    assert edge.cigar == "15643M"
    assert edge.overlap_len == 15643
    assert edge.l1 == 20000
    assert edge.l2 == 30000

    # Dual property access
    dual_prop = edge.dual
    assert dual_prop.from_port == p2.invert()
    assert dual_prop.to_port == p1.invert()
    assert dual_prop.cigar == "15643M"
    assert dual_prop.overlap_len == 15643
    assert dual_prop.l1 == 30000
    assert dual_prop.l2 == 20000

    # Dual method access edge.dual()
    dual_method = edge.dual()
    assert dual_method == dual_prop

    # Invert operator
    assert ~edge == dual_prop

    # Involution invariant: dual(dual(e)) == e
    assert edge.dual.dual == edge
    assert edge.dual().dual() == edge
    assert ~(~edge) == edge


def test_dovetail_edge_immutability_and_hashing():
    """Verify DovetailEdge is frozen and can be indexed in sets/dicts."""
    p1 = SegmentPort(name="utg001", strand=Strand.FORWARD)
    p2 = SegmentPort(name="utg002", strand=Strand.FORWARD)

    e1 = DovetailEdge(p1, p2, "100M", 100)
    e2 = DovetailEdge(p1, p2, "100M", 100)
    assert e1 == e2
    assert hash(e1) == hash(e2)

    edge_set = {e1, e2}
    assert len(edge_set) == 1

    with pytest.raises(FrozenInstanceError):
        e1.overlap_len = 200  # type: ignore[misc]
