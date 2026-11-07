"""rDNA cassette graph-walking resolver with head-to-tail copy chaining.

Implements BIO-001 (unconditional terminal completion) and BIO-005 (lossless dovetail
overlap assembly) for chromosome ends entering nucleolar organizer regions (NOR) and
ribosomal RNA gene clusters. Traverses connected rDNA-bearing unitigs, ensuring
strictly consistent 18S -> 5.8S -> 28S transcriptional orientation without inversion switches,
preventing infinite loops/cycles, chaining maximum consistent rDNA copies, and terminating
at verified distal telomeric caps.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Sequence
from dataclasses import dataclass
from io import StringIO
from pathlib import Path
from typing import TextIO

from gfapdefinisher.annotations import NodeAnnotations
from gfapdefinisher.cigar import cigar_query_consumption
from gfapdefinisher.graph import BidirectedGraph
from gfapdefinisher.ledger import EndEvidence
from gfapdefinisher.models import DovetailEdge, SegmentPort, Strand
from gfapdefinisher.policy import DEFAULT_POLICY, Policy
from gfapdefinisher.quirks import (
    EndState,
    EndStateRecord,
    IntervalQuirk,
    QuirkRecord,
    build_feature_id,
    emit_end_states_tsv,
    emit_quirks_tsv,
)
from gfapdefinisher.resolver import reverse_complement


@dataclass(slots=True)
class RdnaResolution:
    """Resolution result for an rDNA cassette and its distal telomeric cap."""

    chrom: str
    side: str
    state: EndState
    end_state_record: EndStateRecord
    quirk_record: QuirkRecord
    path: list[SegmentPort]
    sequence: str | None = None

    def __iter__(self):
        """Enable tuple unpacking into (end_state_record, quirk_record)."""
        return iter((self.end_state_record, self.quirk_record))

    def emit_end_states_tsv(self, out_file: TextIO | Path | str | None = None) -> str:
        """Emit EndStateRecord to TSV stream or string."""
        return emit_end_states_tsv([self.end_state_record], out_file)

    def emit_quirks_tsv(self, out_file: TextIO | Path | str | None = None) -> str:
        """Emit QuirkRecord to TSV stream or string."""
        buf = StringIO()
        emit_quirks_tsv([self.quirk_record], buf)
        content = buf.getvalue()
        if out_file is not None:
            if isinstance(out_file, str | Path):
                with open(out_file, "w", encoding="utf-8") as f:
                    f.write(content)
            else:
                out_file.write(content)
        return content


class RdnaCassetteResolver:
    """Graph-walking resolver for multi-copy rDNA cassettes and distal telomere caps.

    Per BIO-001, chromosome arms end in authentic telomeric caps, not naked rDNA.
    This resolver traverses collapsed or uncollapsed rDNA tangles, chaining unitigs
    in consistent head-to-tail orientation until reaching the distal telomeric cap.
    """

    def __init__(
        self,
        graph: BidirectedGraph,
        ledger: dict[str, EndEvidence],
        policy: Policy = DEFAULT_POLICY,
        annotations: dict[str, NodeAnnotations] | None = None,
    ) -> None:
        self.graph = graph
        self.ledger = ledger
        self.policy = policy
        self.annotations = annotations

        self._ledger_by_name: dict[str, list[EndEvidence]] = defaultdict(list)
        for ev in self.ledger.values():
            if isinstance(ev, EndEvidence):
                self._ledger_by_name[ev.segment_name].append(ev)

    def _lookup_evidence(self, port: SegmentPort) -> EndEvidence | None:
        """Look up EndEvidence for a port across common indexing keys."""
        if port.to_symbol() in self.ledger:
            return self.ledger[port.to_symbol()]
        end_str = "5p" if port.strand == Strand.MINUS else "3p"
        colon_key = f"{port.name}:{end_str}"
        if colon_key in self.ledger:
            return self.ledger[colon_key]
        strand_key = f"{port.name}:{port.strand.to_symbol()}"
        if strand_key in self.ledger:
            return self.ledger[strand_key]
        if port.name in self.ledger:
            return self.ledger[port.name]
        return None

    def is_rdna_unitig(self, seg_or_port: SegmentPort | str) -> bool:
        """Check whether a segment or port carries ribosomal RNA gene features."""
        name = seg_or_port.name if isinstance(seg_or_port, SegmentPort) else str(seg_or_port).rstrip("+-")

        # Core chromosome arms with single-copy ortholog genes are stems, not repeat unitigs
        if self.annotations is not None and name in self.annotations:
            if len(self.annotations[name].compleasm_hits) > 0:
                return False

        for ev in self._ledger_by_name.get(name, []):
            if ev.usco_hit_count > 0:
                return False

        if self.annotations is not None and name in self.annotations:
            ann = self.annotations[name]
            if bool(ann.barrnap_hits):
                return True

        for ev in self._ledger_by_name.get(name, []):
            if ev.rdna_hit_count > 0 or ev.has_terminal_rdna or bool(ev.rdna_features):
                return True

        return False

    def is_distal_telomere_cap(self, seg_or_port: SegmentPort | str) -> bool:
        """Check whether a segment or port qualifies as a distal telomeric cap."""
        if isinstance(seg_or_port, SegmentPort):
            port = seg_or_port
        else:
            port = self.graph._resolve_port(seg_or_port)

        ev = self._lookup_evidence(port)
        if ev is not None and ev.has_terminal_telomere:
            return True

        for e in self._ledger_by_name.get(port.name, []):
            if e.has_terminal_telomere:
                return True

        return False

    def get_intrinsic_rdna_orientation(self, seg_name: str) -> Strand | None:
        """Determine the 18S -> 5.8S -> 28S transcriptional strand of an rDNA unitig."""
        if self.annotations is not None and seg_name in self.annotations:
            hits = self.annotations[seg_name].barrnap_hits
            h18 = [h for h in hits if "18S" in h.feature_name]
            h28 = [h for h in hits if "28S" in h.feature_name]
            if len(h18) > 0 and len(h28) > 0:
                min_18 = min(h.start for h in h18)
                min_28 = min(h.start for h in h28)
                if min_18 < min_28:
                    return Strand.FORWARD
                elif min_28 < min_18:
                    return Strand.REVERSE

            pos_cnt = sum(1 for h in hits if h.strand == "+")
            neg_cnt = sum(1 for h in hits if h.strand == "-")
            if pos_cnt > neg_cnt:
                return Strand.FORWARD
            elif neg_cnt > pos_cnt:
                return Strand.REVERSE

        for ev in self._ledger_by_name.get(seg_name, []):
            feats = ev.rdna_features
            has_18 = any("18S" in f for f in feats)
            has_28 = any("28S" in f for f in feats)
            if has_18 and not has_28:
                if ev.end == "5p":
                    return Strand.FORWARD
                elif ev.end == "3p":
                    return Strand.REVERSE
            elif has_28 and not has_18:
                if ev.end == "3p":
                    return Strand.FORWARD
                elif ev.end == "5p":
                    return Strand.REVERSE

        return None

    def is_consistent_rdna_orientation(self, port: SegmentPort, previous_ports: Sequence[SegmentPort] = ()) -> bool:
        """Verify that port traversal maintains consistent head-to-tail orientation."""
        intrinsic = self.get_intrinsic_rdna_orientation(port.name)
        if intrinsic is not None and port.strand != intrinsic:
            return False

        return True

    def find_rdna_entry_ports(self) -> list[SegmentPort]:
        """Identify chromosome core stems connecting to rDNA-bearing unitigs."""
        entry_ports: list[SegmentPort] = []
        for seg_name in self.graph.segments:
            if self.is_rdna_unitig(seg_name):
                continue
            for strand in (Strand.FORWARD, Strand.REVERSE):
                port = SegmentPort(seg_name, strand)
                outgoing = self.graph.get_outgoing(port)
                if any(self.is_rdna_unitig(e.to_port) for e in outgoing):
                    entry_ports.append(port)

        entry_ports.sort(key=lambda p: (p.name, p.strand.to_symbol()))
        return entry_ports

    def identify_entry_anchor(self) -> SegmentPort | None:
        """Find the core chromosome arm port connecting into the rDNA tangle."""
        ports = self.find_rdna_entry_ports()
        return ports[0] if len(ports) > 0 else None

    def _find_paths(self, entry_port: SegmentPort) -> tuple[list[list[SegmentPort]], list[list[SegmentPort]]]:
        """Search for valid paths from entry_port through rDNA unitigs to distal caps."""
        max_hops = len(self.graph.segments)
        cap_paths: list[list[SegmentPort]] = []
        rdna_paths: list[list[SegmentPort]] = []

        stack: list[tuple[SegmentPort, list[SegmentPort], set[str]]] = [(entry_port, [entry_port], {entry_port.name})]

        while len(stack) > 0:
            curr_port, path, visited = stack.pop()

            if len(path) > max_hops:
                continue

            outgoing = self.graph.get_outgoing(curr_port)
            advanced = False

            for edge in outgoing:
                nxt = edge.to_port
                if nxt.name in visited:
                    continue

                if self.is_distal_telomere_cap(nxt):
                    new_path = list(path)
                    new_path.append(nxt)
                    cap_paths.append(new_path)
                    advanced = True
                    continue

                if self.is_rdna_unitig(nxt):
                    if not self.is_consistent_rdna_orientation(nxt, path):
                        continue

                    new_path = list(path)
                    new_path.append(nxt)
                    new_visited = set(visited)
                    new_visited.add(nxt.name)
                    stack.append((nxt, new_path, new_visited))
                    advanced = True

            if not advanced and len(path) > 1:
                rdna_paths.append(path)

        return cap_paths, rdna_paths

    def reach_distal_telomere_cap(self, entry_port: SegmentPort | str | None = None) -> list[SegmentPort] | None:
        """Identify the path that terminates at a distal telomeric cap unitig."""
        resolved_entry: SegmentPort | None
        if entry_port is None:
            entries = self.find_rdna_entry_ports()
            resolved_entry = entries[0] if len(entries) > 0 else None
        else:
            resolved_entry = entry_port if isinstance(entry_port, SegmentPort) else self.graph._resolve_port(entry_port)

        if resolved_entry is None:
            return None

        cap_paths, _ = self._find_paths(resolved_entry)
        if len(cap_paths) == 0:
            return None

        def score_path(p: list[SegmentPort]) -> tuple[int, int]:
            rdna_count = sum(1 for port in p if self.is_rdna_unitig(port))
            return (rdna_count, self.compute_path_length(p))

        cap_paths.sort(key=score_path, reverse=True)
        return cap_paths[0]

    def connect_to_telomere_cap(self, entry_port: SegmentPort | str | None = None) -> list[SegmentPort] | None:
        """Alias for reach_distal_telomere_cap."""
        return self.reach_distal_telomere_cap(entry_port)

    def traverse_rdna_cassette(
        self,
        entry_port: SegmentPort | str | None = None,
        include_anchor: bool = False,
    ) -> list[SegmentPort]:
        """Traverse through connected rDNA unitigs chaining consistent copies."""
        resolved_entry: SegmentPort | None
        if entry_port is None:
            entries = self.find_rdna_entry_ports()
            resolved_entry = entries[0] if len(entries) > 0 else None
        else:
            resolved_entry = entry_port if isinstance(entry_port, SegmentPort) else self.graph._resolve_port(entry_port)

        if resolved_entry is None:
            return []

        cap_paths, rdna_paths = self._find_paths(resolved_entry)

        best_path: list[SegmentPort] = []
        if len(cap_paths) > 0:

            def score_cap(p: list[SegmentPort]) -> tuple[int, int]:
                cnt = sum(1 for port in p if self.is_rdna_unitig(port))
                return (cnt, self.compute_path_length(p))

            cap_paths.sort(key=score_cap, reverse=True)
            chosen = cap_paths[0]
            if include_anchor:
                best_path = chosen
            else:
                best_path = [p for p in chosen if self.is_rdna_unitig(p)]
        elif len(rdna_paths) > 0:

            def score_rdna(p: list[SegmentPort]) -> tuple[int, int]:
                cnt = sum(1 for port in p if self.is_rdna_unitig(port))
                return (cnt, self.compute_path_length(p))

            rdna_paths.sort(key=score_rdna, reverse=True)
            chosen = rdna_paths[0]
            if include_anchor:
                best_path = chosen
            else:
                best_path = [p for p in chosen if self.is_rdna_unitig(p)]

        return best_path

    def traverse_head_to_tail(
        self,
        entry_port: SegmentPort | str | None = None,
        include_anchor: bool = False,
    ) -> list[SegmentPort]:
        """Alias for traverse_rdna_cassette."""
        return self.traverse_rdna_cassette(entry_port, include_anchor=include_anchor)

    def compute_path_length(self, path: Sequence[SegmentPort]) -> int:
        """Compute exact dovetail assembled length of a path in base pairs."""
        if len(path) == 0:
            return 0

        total_bp = 0
        for idx, port in enumerate(path):
            seg = self.graph.segments.get(port.name)
            seg_len = seg.length if seg is not None else 0
            if idx == 0:
                total_bp += seg_len
            else:
                prev_port = path[idx - 1]
                overlap = 0
                for edge in self.graph.get_outgoing(prev_port):
                    if edge.to_port == port:
                        overlap = edge.overlap_len
                        break
                total_bp += max(0, seg_len - overlap)

        return total_bp

    def assemble_contiguous_path(self, path: Sequence[SegmentPort] | None = None) -> str | None:
        """Stitch traversed unitigs into a contiguous sequence using exact dovetails (BIO-005)."""
        target_path = path
        if target_path is None:
            target_path = self.reach_distal_telomere_cap()
        if target_path is None or len(target_path) == 0:
            return None

        first_port = target_path[0]
        first_seg = self.graph.segments.get(first_port.name)
        if first_seg is None or first_seg.sequence is None:
            return None

        seq = first_seg.sequence if first_port.strand == Strand.FORWARD else reverse_complement(first_seg.sequence)

        assembled_chunks: list[str] = [seq]

        for idx in range(1, len(target_path)):
            prev_port = target_path[idx - 1]
            curr_port = target_path[idx]

            curr_seg = self.graph.segments.get(curr_port.name)
            if curr_seg is None or curr_seg.sequence is None:
                return None

            curr_seq = (
                curr_seg.sequence if curr_port.strand == Strand.FORWARD else reverse_complement(curr_seg.sequence)
            )

            edge: DovetailEdge | None = None
            for e in self.graph.get_outgoing(prev_port):
                if e.to_port == curr_port:
                    edge = e
                    break

            if edge is not None and edge.cigar and edge.cigar != "*":
                q_consume = cigar_query_consumption(edge.cigar)
            elif edge is not None:
                q_consume = edge.overlap_len
            else:
                q_consume = 0

            assembled_chunks.append(curr_seq[q_consume:])

        return "".join(assembled_chunks)

    def resolve_rdna_cassette(
        self,
        chrom: str = "II",
        side: str = "5p",
        entry_port: SegmentPort | str | None = None,
        gap_index: int = 1,
        notes: str = "",
    ) -> RdnaResolution:
        """Assemble traversed rDNA path and emit EndStateRecord and QuirkRecord."""
        resolved_entry: SegmentPort | None
        if entry_port is None:
            entries = self.find_rdna_entry_ports()
            resolved_entry = entries[0] if len(entries) > 0 else None
        else:
            resolved_entry = entry_port if isinstance(entry_port, SegmentPort) else self.graph._resolve_port(entry_port)

        path_to_cap = self.reach_distal_telomere_cap(resolved_entry) if resolved_entry is not None else None

        if path_to_cap is not None and len(path_to_cap) > 1:
            cap_port = path_to_cap[-1]
            assembled_seq = self.assemble_contiguous_path(path_to_cap)
            path_len = self.compute_path_length(path_to_cap)
            all_unitigs = tuple(p.name for p in path_to_cap)
            rdna_unitigs = tuple(p.name for p in path_to_cap if self.is_rdna_unitig(p))

            feature_id_end = build_feature_id(EndState.END_SHARED_EXEMPT, chrom, side)
            feature_id_quirk = build_feature_id(IntervalQuirk.RDNA_CASSETTE, chrom, side, index=gap_index)

            quirk_notes = (
                f"rDNA cassette chained {len(rdna_unitigs)} copies terminating at distal cap {cap_port.to_symbol()}"
            )
            if notes:
                quirk_notes = f"{quirk_notes}; {notes}"

            end_notes = f"Shared exempt rDNA cassette resolved to distal telomeric cap {cap_port.to_symbol()}"
            if notes:
                end_notes = f"{end_notes}; {notes}"

            core_cov: float | None = None
            anchor_seg = self.graph.segments.get(path_to_cap[0].name)
            if anchor_seg is not None and anchor_seg.depth > 0:
                core_cov = anchor_seg.depth

            quirk_record = QuirkRecord(
                feature_id=feature_id_quirk,
                chrom=chrom,
                side=side,
                quirk_type=IntervalQuirk.RDNA_CASSETTE,
                length_min_bp=path_len,
                length_median_bp=path_len,
                length_max_bp=path_len,
                left_anchor_port=path_to_cap[0].to_symbol(),
                right_anchor_port=cap_port.to_symbol(),
                agp_component_type="W",
                agp_linkage_evidence="graph-walk",
                unitig_ids=all_unitigs,
                core_coverage_bp_median=core_cov,
                representative_sequence_selected=assembled_seq is not None,
                notes=quirk_notes,
            )

            end_state_record = EndStateRecord(
                feature_id=feature_id_end,
                chrom=chrom,
                side=side,
                state=EndState.END_SHARED_EXEMPT,
                somatic_port=None,
                germline_port=None,
                shared_port=cap_port.to_symbol(),
                core_coverage_bp_median=core_cov,
                resolved=True,
                notes=end_notes,
            )

            return RdnaResolution(
                chrom=chrom,
                side=side,
                state=EndState.END_SHARED_EXEMPT,
                end_state_record=end_state_record,
                quirk_record=quirk_record,
                path=path_to_cap,
                sequence=assembled_seq,
            )

        # Distal cap not reached
        feature_id_end = build_feature_id(EndState.END_SHARED_EXEMPT, chrom, side)
        feature_id_quirk = build_feature_id(IntervalQuirk.RDNA_CASSETTE, chrom, side, index=gap_index)

        unresolved_path = (
            self.traverse_rdna_cassette(resolved_entry, include_anchor=True) if resolved_entry is not None else []
        )
        all_unitigs = tuple(p.name for p in unresolved_path)
        path_len = self.compute_path_length(unresolved_path) if len(unresolved_path) > 0 else None
        assembled_seq = self.assemble_contiguous_path(unresolved_path) if len(unresolved_path) > 0 else None

        left_str = resolved_entry.to_symbol() if resolved_entry is not None else None

        unresolved_quirk_notes = "Unresolved rDNA cassette tangle: distal telomere cap not reached"
        if notes:
            unresolved_quirk_notes = f"{unresolved_quirk_notes}; {notes}"

        unresolved_end_notes = "Unresolved rDNA cassette end"
        if notes:
            unresolved_end_notes = f"{unresolved_end_notes}; {notes}"

        quirk_record = QuirkRecord(
            feature_id=feature_id_quirk,
            chrom=chrom,
            side=side,
            quirk_type=IntervalQuirk.RDNA_CASSETTE,
            length_min_bp=path_len,
            length_median_bp=path_len,
            length_max_bp=path_len,
            left_anchor_port=left_str,
            right_anchor_port=None,
            agp_component_type="U",
            agp_linkage_evidence="unresolved-tangle",
            unitig_ids=all_unitigs,
            representative_sequence_selected=assembled_seq is not None,
            notes=unresolved_quirk_notes,
        )

        end_state_record = EndStateRecord(
            feature_id=feature_id_end,
            chrom=chrom,
            side=side,
            state=EndState.END_SHARED_EXEMPT,
            shared_port=None,
            resolved=False,
            notes=unresolved_end_notes,
        )

        return RdnaResolution(
            chrom=chrom,
            side=side,
            state=EndState.END_SHARED_EXEMPT,
            end_state_record=end_state_record,
            quirk_record=quirk_record,
            path=unresolved_path,
            sequence=assembled_seq,
        )
