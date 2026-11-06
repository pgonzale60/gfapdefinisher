"""Lossless GFA1 parser and bidirected assembly graph core.

Provides bidirected graph representations with Watson-Crick duality invariance,
lossless GFA1 parsing preserving all tags and A lines, and degree/tip queries.
"""

from __future__ import annotations

import gzip
from pathlib import Path
from typing import Any

from gfapdefinisher.cigar import cigar_overlap_length
from gfapdefinisher.models import DovetailEdge, Segment, SegmentPort, Strand


def _parse_gfa_tag(tag_str: str) -> tuple[str, Any]:
    """Parse a GFA1 tag string formatted as TAG:TYPE:VALUE.

    Supports integer (:i:), float (:f:), and string/char types losslessly.
    """
    try:
        tag_name, tag_type, tag_val = tag_str.split(":", 2)
    except ValueError:
        return tag_str, tag_str

    if tag_type == "i":
        try:
            return tag_name, int(tag_val)
        except ValueError:
            return tag_name, tag_val
    elif tag_type == "f":
        try:
            return tag_name, float(tag_val)
        except ValueError:
            return tag_name, tag_val
    else:
        return tag_name, tag_val


class BidirectedGraph:
    """Bidirected graph core for assembly unitigs and dovetail overlaps.

    Maintains segments and strand-aware dovetail edges with strict Watson-Crick
    duality invariance: every edge u(strand) -> v(strand) has an automatic dual
    ~v(strand) -> ~u(strand). Preserves all GFA tags and alignment records losslessly.
    """

    def __init__(self) -> None:
        self.segments: dict[str, Segment] = {}
        self.edges: dict[SegmentPort, list[DovetailEdge]] = {}
        self._incoming: dict[SegmentPort, list[DovetailEdge]] = {}
        self.segment_tags: dict[str, dict[str, Any]] = {}
        self.raw_segment_tags: dict[str, list[str]] = {}
        self.edge_tags: dict[tuple[SegmentPort, SegmentPort], dict[str, Any]] = {}
        self.raw_edge_tags: dict[tuple[SegmentPort, SegmentPort], list[str]] = {}
        self.a_lines: list[str] = []
        self.headers: list[str] = []
        self.paths: list[str] = []
        self.other_lines: list[str] = []
        self.gfa_path: Path | None = None

    @property
    def tags(self) -> dict[str, dict[str, Any]]:
        """Access segment tags mapping."""
        return self.segment_tags

    @staticmethod
    def _resolve_port(port: SegmentPort | str, strand: Strand | str | None = None) -> SegmentPort:
        """Resolve a port representation to a canonical SegmentPort."""
        if isinstance(port, SegmentPort):
            if strand is not None:
                return SegmentPort(name=port.name, strand=strand)
            return port
        if strand is not None:
            return SegmentPort(name=port, strand=strand)
        if port.endswith("+"):
            return SegmentPort(name=port[:-1], strand=Strand.FORWARD)
        if port.endswith("-"):
            return SegmentPort(name=port[:-1], strand=Strand.REVERSE)
        return SegmentPort(name=port, strand=Strand.FORWARD)

    def add_segment(
        self,
        name: str,
        length: int,
        depth: float = 0.0,
        sequence: str | None = None,
        tags: dict[str, Any] | None = None,
        raw_tags: list[str] | None = None,
    ) -> Segment:
        """Register a segment with length, depth, and optional sequence and tags."""
        seg = Segment(name=name, length=length, depth=depth, sequence=sequence)
        self.segments[name] = seg

        # Pre-initialize port edge buckets for both strands
        fwd_port = SegmentPort(name=name, strand=Strand.FORWARD)
        rev_port = SegmentPort(name=name, strand=Strand.REVERSE)
        self.edges.setdefault(fwd_port, [])
        self.edges.setdefault(rev_port, [])
        self._incoming.setdefault(fwd_port, [])
        self._incoming.setdefault(rev_port, [])

        if tags is not None:
            self.segment_tags[name] = tags
        elif name not in self.segment_tags:
            self.segment_tags[name] = {}

        if raw_tags is not None:
            self.raw_segment_tags[name] = raw_tags
        elif name not in self.raw_segment_tags:
            self.raw_segment_tags[name] = []

        return seg

    def _record_directed_edge(
        self,
        edge: DovetailEdge,
        tags: dict[str, Any] | None = None,
        raw_tags: list[str] | None = None,
    ) -> None:
        """Record a single directed edge in outgoing and incoming registries."""
        port_edges = self.edges.setdefault(edge.from_port, [])
        if edge not in port_edges:
            port_edges.append(edge)

        inc_edges = self._incoming.setdefault(edge.to_port, [])
        if edge not in inc_edges:
            inc_edges.append(edge)

        key = (edge.from_port, edge.to_port)
        if tags is not None:
            self.edge_tags[key] = tags
        elif key not in self.edge_tags:
            self.edge_tags[key] = {}

        if raw_tags is not None:
            self.raw_edge_tags[key] = raw_tags
        elif key not in self.raw_edge_tags:
            self.raw_edge_tags[key] = []

    def add_edge(
        self,
        from_port: SegmentPort | str,
        to_port: SegmentPort | str,
        cigar: str,
        l1: int = 0,
        l2: int = 0,
        add_dual: bool = True,
        tags: dict[str, Any] | None = None,
        raw_tags: list[str] | None = None,
    ) -> DovetailEdge:
        """Add a dovetail overlap edge and automatically maintain its reverse dual."""
        p_from = self._resolve_port(from_port)
        p_to = self._resolve_port(to_port)

        overlap_len = cigar_overlap_length(cigar)
        if l1 == 0 and p_from.name in self.segments:
            l1 = self.segments[p_from.name].length
        if l2 == 0 and p_to.name in self.segments:
            l2 = self.segments[p_to.name].length

        edge = DovetailEdge(
            from_port=p_from,
            to_port=p_to,
            cigar=cigar,
            overlap_len=overlap_len,
            l1=l1,
            l2=l2,
        )

        self._record_directed_edge(edge, tags, raw_tags)

        if add_dual:
            dual_edge = edge.dual
            self._record_directed_edge(dual_edge, tags, raw_tags)

        return edge

    def get_outgoing(self, port: SegmentPort | str, strand: Strand | str | None = None) -> list[DovetailEdge]:
        """Return list of outgoing dovetail edges extending from port."""
        resolved = self._resolve_port(port, strand)
        return list(self.edges.get(resolved, []))

    def get_incoming(self, port: SegmentPort | str, strand: Strand | str | None = None) -> list[DovetailEdge]:
        """Return list of incoming dovetail edges entering port."""
        resolved = self._resolve_port(port, strand)
        return list(self._incoming.get(resolved, []))

    def degree(self, port: SegmentPort | str, strand: Strand | str | None = None) -> int:
        """Return the outgoing degree from the given oriented port."""
        resolved = self._resolve_port(port, strand)
        return len(self.edges.get(resolved, []))

    def is_tip(self, port: SegmentPort | str, strand: Strand | str | None = None) -> bool:
        """Return True if port is a tip (has outgoing degree of 0)."""
        return self.degree(port, strand) == 0

    @property
    def canonical_edges(self) -> list[DovetailEdge]:
        """Return deduplicated canonical dovetail edges (one per dual pair)."""
        seen: set[tuple[SegmentPort, SegmentPort]] = set()
        canonical: list[DovetailEdge] = []
        for edge_list in self.edges.values():
            for edge in edge_list:
                pair = (edge.from_port, edge.to_port)
                dual_pair = (edge.to_port.invert(), edge.from_port.invert())
                if pair not in seen and dual_pair not in seen:
                    seen.add(pair)
                    canonical.append(edge)
        return canonical

    def load_all_sequences(self, gfa_path: str | Path | None = None) -> None:
        """Load sequences from GFA file into segments that lack sequence in memory."""
        target_path = Path(gfa_path) if gfa_path is not None else self.gfa_path
        if target_path is None:
            raise ValueError("No GFA path available to load sequences from.")

        is_gz = target_path.name.endswith(".gz")
        opener = gzip.open if is_gz else open
        with opener(target_path, "rt", encoding="utf-8") as f:
            for line in f:
                if not line.startswith("S\t"):
                    continue
                parts = line.rstrip("\r\n").split("\t")
                try:
                    _rec, seg_name, raw_seq, *_ = parts
                except ValueError:
                    continue
                if raw_seq != "*" and seg_name in self.segments:
                    self.segments[seg_name].sequence = raw_seq

    def get_sequence(self, name: str) -> str | None:
        """Get sequence for segment, loading on-demand if missing and gfa_path set."""
        if name not in self.segments:
            return None
        seg = self.segments[name]
        if seg.sequence is not None:
            return seg.sequence
        if self.gfa_path is not None:
            self.load_all_sequences()
            return seg.sequence
        return None

    def write_gfa(self, out_path: str | Path) -> None:
        """Write graph to GFA1 losslessly preserving segments, links, tags, and A lines."""
        path = Path(out_path)
        is_gz = path.name.endswith(".gz")
        opener = gzip.open if is_gz else open
        with opener(path, "wt", encoding="utf-8") as f:
            for header in self.headers:
                f.write(f"{header}\n")
            if not self.headers:
                f.write("H\tVN:Z:1.0\n")

            for name, seg in self.segments.items():
                seq_str = seg.sequence if seg.sequence is not None else "*"
                tag_strs = list(self.raw_segment_tags.get(name, []))
                if not tag_strs:
                    tag_strs.append(f"LN:i:{seg.length}")
                    if seg.depth != 0.0:
                        tag_strs.append(f"rd:i:{int(seg.depth)}")
                tag_part = "\t" + "\t".join(tag_strs) if tag_strs else ""
                f.write(f"S\t{name}\t{seq_str}{tag_part}\n")

            written_pairs: set[tuple[SegmentPort, SegmentPort]] = set()
            for edge_list in self.edges.values():
                for edge in edge_list:
                    pair = (edge.from_port, edge.to_port)
                    dual_pair = (edge.to_port.invert(), edge.from_port.invert())
                    if pair in written_pairs or dual_pair in written_pairs:
                        continue
                    written_pairs.add(pair)
                    tag_strs = list(self.raw_edge_tags.get(pair, []))
                    if not tag_strs:
                        if edge.l1 != 0:
                            tag_strs.append(f"L1:i:{edge.l1}")
                        if edge.l2 != 0:
                            tag_strs.append(f"L2:i:{edge.l2}")
                    tag_part = "\t" + "\t".join(tag_strs) if tag_strs else ""
                    f.write(
                        f"L\t{edge.from_port.name}\t{edge.from_port.strand.to_symbol()}\t"
                        f"{edge.to_port.name}\t{edge.to_port.strand.to_symbol()}\t"
                        f"{edge.cigar}{tag_part}\n"
                    )

            for a_line in self.a_lines:
                f.write(f"{a_line}\n")

    def __len__(self) -> int:
        return len(self.segments)

    def __contains__(self, name: str) -> bool:
        return name in self.segments

    def __getitem__(self, name: str) -> Segment:
        return self.segments[name]


def load_gfa(
    gfa_path: str | Path,
    load_sequences: bool = False,
    a_lines_store: list[str] | None = None,
) -> BidirectedGraph:
    """Parse a GFA1 file into a BidirectedGraph losslessly.

    Preserves segment lengths (LN:i:), depths (rd:i: or dp:f:), links,
    tags, and alignment (A) lines. Supports optional sequence loading.

    Args:
        gfa_path: Path to GFA1 file (supports .gfa and .gfa.gz).
        load_sequences: If True, reads nucleotide sequences into memory.
        a_lines_store: Optional list to append raw GFA A lines to.

    Returns:
        Populated BidirectedGraph.
    """
    path = Path(gfa_path)
    graph = BidirectedGraph()
    graph.gfa_path = path

    is_gz = path.name.endswith(".gz")
    opener = gzip.open if is_gz else open

    with opener(path, "rt", encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            record_type = parts[0]

            if record_type == "H":
                graph.headers.append(line)

            elif record_type == "S":
                try:
                    _rec, seg_name, raw_seq, *raw_tag_list = parts
                except ValueError:
                    continue

                tags: dict[str, Any] = {}
                for tag_str in raw_tag_list:
                    tag_name, tag_val = _parse_gfa_tag(tag_str)
                    tags[tag_name] = tag_val

                length = 0
                if "LN" in tags and isinstance(tags["LN"], int):
                    length = tags["LN"]
                elif raw_seq != "*":
                    length = len(raw_seq)

                depth = 0.0
                if "rd" in tags:
                    try:
                        depth = float(tags["rd"])
                    except (ValueError, TypeError):
                        pass
                elif "dp" in tags:
                    try:
                        depth = float(tags["dp"])
                    except (ValueError, TypeError):
                        pass

                sequence = None
                if load_sequences and raw_seq != "*":
                    sequence = raw_seq

                graph.add_segment(
                    name=seg_name,
                    length=length,
                    depth=depth,
                    sequence=sequence,
                    tags=tags,
                    raw_tags=raw_tag_list,
                )

            elif record_type == "L":
                try:
                    _rec, from_name, from_strand_sym, to_name, to_strand_sym, cigar, *raw_tag_list = parts
                except ValueError:
                    continue

                from_strand = Strand.from_symbol(from_strand_sym)
                to_strand = Strand.from_symbol(to_strand_sym)
                from_port = SegmentPort(name=from_name, strand=from_strand)
                to_port = SegmentPort(name=to_name, strand=to_strand)

                tags = {}
                for tag_str in raw_tag_list:
                    tag_name, tag_val = _parse_gfa_tag(tag_str)
                    tags[tag_name] = tag_val

                l1 = tags.get("L1", 0)
                if not isinstance(l1, int):
                    l1 = 0
                l2 = tags.get("L2", 0)
                if not isinstance(l2, int):
                    l2 = 0

                graph.add_edge(
                    from_port=from_port,
                    to_port=to_port,
                    cigar=cigar,
                    l1=l1,
                    l2=l2,
                    add_dual=True,
                    tags=tags,
                    raw_tags=raw_tag_list,
                )

            elif record_type == "A":
                graph.a_lines.append(line)
                if a_lines_store is not None:
                    a_lines_store.append(line)

            elif record_type == "P":
                graph.paths.append(line)

            else:
                graph.other_lines.append(line)

    return graph
