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
    distal_cap_port: SegmentPort | None = None
    capped_via_heuristic: bool = False
    reached_distal_cap: bool = False

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
        self.ledger_by_port: dict[str, EndEvidence] = ledger
        self.policy = policy
        self.annotations = annotations

        self._ledger_by_name: dict[str, list[EndEvidence]] = defaultdict(list)
        for ev in self.ledger.values():
            if isinstance(ev, EndEvidence):
                self._ledger_by_name[ev.segment_name].append(ev)

    def _is_contaminant(self, seg_name: str) -> bool:
        """Check whether a segment is flagged as contaminant or non-eukaryotic."""
        clean_name = seg_name.rstrip("+-") if isinstance(seg_name, str) else seg_name.name
        contam_classes = {
            "bacteria",
            "archaea",
            "prokarya",
            "contaminant",
            "organelle",
            "mitochondrion",
            "plastid",
        }
        if self.annotations is not None and clean_name in self.annotations:
            ann = self.annotations[clean_name]
            if getattr(ann, "is_contaminant", None) and ann.is_contaminant():
                return True
            tiara_cls = getattr(ann, "tiara_classification", None)
            if tiara_cls and str(tiara_cls).lower().strip() in contam_classes:
                return True

        port_keys = (
            f"{clean_name}+",
            f"{clean_name}-",
            f"{clean_name}:5p",
            f"{clean_name}:3p",
            clean_name,
        )
        for key in port_keys:
            ev = self.ledger_by_port.get(key)
            if ev is not None:
                if getattr(ev, "is_contaminant", False):
                    return True
                cls_str = str(getattr(ev, "tiara_class", "")).lower().strip()
                if cls_str in contam_classes:
                    return True

        for ev in self._ledger_by_name.get(clean_name, []):
            if getattr(ev, "is_contaminant", False):
                return True
            cls_str = str(getattr(ev, "tiara_class", "")).lower().strip()
            if cls_str in contam_classes:
                return True

        return False

    def _is_contaminant_unitig(self, seg_name: str) -> bool:
        """Check whether a unitig is flagged as contaminant (alias for _is_contaminant)."""
        return self._is_contaminant(seg_name)

    def _get_usco_count(self, seg_name: str) -> int:
        """Retrieve total USCO hit count for segment."""
        clean_name = seg_name.rstrip("+-") if isinstance(seg_name, str) else seg_name.name
        count = 0
        if self.annotations is not None:
            if hasattr(self.annotations, "usco_counts"):
                count += self.annotations.usco_counts.get(clean_name, 0)
            elif hasattr(self.annotations, "get") and clean_name in self.annotations:
                ann = self.annotations[clean_name]
                if hasattr(ann, "compleasm_hits"):
                    count += len(ann.compleasm_hits)
        for ev in self._ledger_by_name.get(clean_name, []):
            if ev.usco_hit_count > count:
                count = ev.usco_hit_count
        return count

    def _has_terminal_telomere(self, port: SegmentPort) -> bool:
        """Check whether an oriented port carries a terminal telomere."""
        ev = self._lookup_evidence(port)
        if ev is not None:
            return bool(ev.has_terminal_telomere)
        end_str = "5p" if port.strand == Strand.MINUS else "3p"
        for e in self._ledger_by_name.get(port.name, []):
            if e.end == end_str:
                return bool(e.has_terminal_telomere)
        return False

    def _get_rdna_hit_count(self, seg_id: str) -> int:
        """Return the number of rDNA hits for a segment."""
        if self.annotations is not None:
            if hasattr(self.annotations, "rdna_hits"):
                hits_map = self.annotations.rdna_hits
                if hasattr(hits_map, "get"):
                    return len(hits_map.get(seg_id, ()))
            if hasattr(self.annotations, "get") and seg_id in self.annotations:
                ann = self.annotations[seg_id]
                if hasattr(ann, "barrnap_hits"):
                    return len(ann.barrnap_hits)
        for ev in self._ledger_by_name.get(seg_id, []):
            if ev.rdna_hit_count > 0:
                return ev.rdna_hit_count
        return 0

    def _has_any_usco(self) -> bool:
        """Check if any segment in the graph or ledger has USCO hits."""
        if self.annotations is not None:
            for ann in self.annotations.values():
                if len(ann.compleasm_hits) > 0:
                    return True
        for ev_list in self._ledger_by_name.values():
            for ev in ev_list:
                if ev.usco_hit_count > 0:
                    return True
        if hasattr(self, "ledger_by_port") and self.ledger_by_port:
            for ev in self.ledger_by_port.values():
                if isinstance(ev, EndEvidence) and ev.usco_hit_count > 0:
                    return True
        return False

    def _get_rdna_segments(self) -> set[str]:
        """Return the set of valid eukaryotic rDNA-bearing segment names."""
        rdna_segs: set[str] = set()
        for seg_name in self.graph.segments:
            if self._is_contaminant(seg_name):
                continue
            if self.is_rdna_unitig(seg_name):
                rdna_segs.add(seg_name)
        return rdna_segs

    def _is_telomeric_cap(self, port: SegmentPort) -> bool:
        """Check whether a port qualifies as a distal telomeric cap."""
        if self._is_contaminant(port.name):
            return False
        if self._get_usco_count(port.name) > 0:
            return False
        ev = self._lookup_evidence(port)
        if ev is not None:
            return bool(ev.has_terminal_telomere)
        return False

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

        if self._is_contaminant(name):
            return False

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
        name = seg_or_port.name if isinstance(seg_or_port, SegmentPort) else str(seg_or_port).rstrip("+-")

        if self._is_contaminant(name):
            return False
        if self._get_usco_count(name) > 0:
            return False

        if isinstance(seg_or_port, SegmentPort):
            port = seg_or_port
            ev = self._lookup_evidence(port)
            if ev is not None:
                return bool(ev.has_terminal_telomere)
            return False

        if seg_or_port.endswith(("+", "-")):
            port = self.graph._resolve_port(seg_or_port)
            return self.is_distal_telomere_cap(port)

        seg_name = seg_or_port
        for e in self._ledger_by_name.get(seg_name, []):
            if e.has_terminal_telomere:
                return True

        return False

    def get_intrinsic_rdna_orientation(self, seg_name: str) -> Strand | None:
        """Determine the 18S -> 5.8S -> 28S transcriptional strand of an rDNA unitig."""
        if self.annotations is not None and seg_name in self.annotations:
            hits = self.annotations[seg_name].barrnap_hits
            pos_hits = [h for h in hits if h.strand == "+"]
            neg_hits = [h for h in hits if h.strand == "-"]

            # Check positive strand 18S -> 28S structure (increasing coordinates)
            pos_18 = [h for h in pos_hits if "18S" in h.feature_name]
            pos_28 = [h for h in pos_hits if "28S" in h.feature_name]
            pos_ordered = False
            pos_inverted = False
            if len(pos_18) > 0 and len(pos_28) > 0:
                min_18_pos = min(h.start for h in pos_18)
                min_28_pos = min(h.start for h in pos_28)
                if min_18_pos < min_28_pos:
                    pos_ordered = True
                elif min_28_pos < min_18_pos:
                    pos_inverted = True

            # Check negative strand 18S -> 28S structure (decreasing coordinates in unitig)
            neg_18 = [h for h in neg_hits if "18S" in h.feature_name]
            neg_28 = [h for h in neg_hits if "28S" in h.feature_name]
            neg_ordered = False
            neg_inverted = False
            if len(neg_18) > 0 and len(neg_28) > 0:
                min_18_neg = min(h.start for h in neg_18)
                min_28_neg = min(h.start for h in neg_28)
                if min_18_neg > min_28_neg:
                    neg_ordered = True
                elif min_28_neg > min_18_neg:
                    neg_inverted = True

            # Internal inversion on either strand cannot represent valid transcriptional tandem unit
            if pos_inverted or neg_inverted:
                return None

            if pos_ordered and not neg_ordered:
                return Strand.FORWARD
            if neg_ordered and not pos_ordered:
                return Strand.REVERSE

            # Majority strand fallback when only one subunit is present
            pos_cnt = len(pos_hits)
            neg_cnt = len(neg_hits)
            if pos_cnt > 0 and neg_cnt == 0:
                return Strand.FORWARD
            elif neg_cnt > 0 and pos_cnt == 0:
                return Strand.REVERSE
            elif pos_cnt > neg_cnt:
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

    def is_consistent_rdna_orientation(
        self,
        port: SegmentPort,
        previous_ports: Sequence[SegmentPort] | Strand | None = (),
    ) -> bool:
        """Verify that port traversal maintains consistent head-to-tail orientation."""
        if not self.is_rdna_unitig(port):
            return False

        intrinsic = self.get_intrinsic_rdna_orientation(port.name)
        if intrinsic is not None:
            if port.strand != intrinsic:
                return False
        else:
            # If unitig carries Barrnap rDNA annotations but no coherent orientation, reject inversion/chimera
            if self.annotations is not None:
                if hasattr(self.annotations, "get") and port.name in self.annotations:
                    ann = self.annotations[port.name]
                    if hasattr(ann, "barrnap_hits") and len(ann.barrnap_hits) > 0:
                        return False
                elif hasattr(self.annotations, "rdna_hits"):
                    hits_map = self.annotations.rdna_hits
                    if hasattr(hits_map, "get") and len(hits_map.get(port.name, ())) > 0:
                        return False

        if isinstance(previous_ports, Strand):
            if port.strand != previous_ports:
                return False
        elif previous_ports and isinstance(previous_ports, Sequence) and not isinstance(previous_ports, str | bytes):
            prev_rdna = [p for p in previous_ports if self.is_rdna_unitig(p)]
            if len(prev_rdna) > 0:
                first_intrinsic = self.get_intrinsic_rdna_orientation(prev_rdna[0].name)
                if first_intrinsic is not None and intrinsic is not None:
                    if intrinsic != first_intrinsic:
                        return False

        return True

    def find_heuristic_rdna_cap(
        self,
        expected_dir: Strand | None,
        visited_ids: set[str],
    ) -> SegmentPort | None:
        """Find an orphan rDNA-bearing telomeric cap unitig matching orientation and purity criteria."""
        candidates: list[SegmentPort] = []
        for seg_id in sorted(self.graph.segments.keys()):
            if seg_id in visited_ids:
                continue
            if self._is_contaminant_unitig(seg_id):
                continue
            if not self.is_rdna_unitig(seg_id):
                continue

            p_fwd = SegmentPort(seg_id, Strand.FORWARD)
            p_rev = SegmentPort(seg_id, Strand.REVERSE)
            if self.graph.degree(p_fwd) != 0 or self.graph.degree(p_rev) != 0:
                continue

            if self._get_usco_count(seg_id) != 0:
                continue

            has_fwd = self._has_terminal_telomere(p_fwd)
            has_rev = self._has_terminal_telomere(p_rev)
            if has_fwd == has_rev:
                continue

            cand_port = SegmentPort(seg_id, Strand.FORWARD if has_fwd else Strand.REVERSE)
            if not self.is_consistent_rdna_orientation(cand_port, expected_dir):
                continue

            candidates.append(cand_port)

        if len(candidates) == 0:
            return None

        candidates.sort(
            key=lambda p: (
                -self._get_rdna_hit_count(p.name),
                -(self.graph.segments[p.name].length if p.name in self.graph.segments else 0),
                p.name,
            )
        )
        return candidates[0]

    def _find_entry_stems(self) -> list[SegmentPort]:
        """Identify chromosome core stems connecting to rDNA-bearing unitigs."""
        rdna_segments = self._get_rdna_segments()
        has_usco = self._has_any_usco()
        entry_ports: list[SegmentPort] = []

        for seg_name in self.graph.segments:
            if self._is_contaminant(seg_name):
                continue
            if seg_name in rdna_segments:
                continue
            if self.is_distal_telomere_cap(seg_name):
                continue
            if has_usco and self._get_usco_count(seg_name) == 0:
                continue

            for strand in (Strand.FORWARD, Strand.REVERSE):
                port = SegmentPort(seg_name, strand)
                outgoing = self.graph.get_outgoing(port)
                if any(e.to_port.name in rdna_segments for e in outgoing):
                    entry_ports.append(port)

        def rank_entry(p: SegmentPort) -> tuple[int, int, str, str]:
            usco_flag = 1 if self._get_usco_count(p.name) > 0 else 0
            seg = self.graph.segments.get(p.name)
            seg_len = seg.length if seg is not None else 0
            # Sort descending by USCO presence and length, with deterministic name/strand
            return (-usco_flag, -seg_len, p.name, p.strand.to_symbol())

        entry_ports.sort(key=rank_entry)
        return entry_ports

    def find_rdna_entry_ports(self) -> list[SegmentPort]:
        """Identify chromosome core stems connecting to rDNA-bearing unitigs."""
        return self._find_entry_stems()

    def identify_entry_anchor(self) -> SegmentPort | None:
        """Find the core chromosome arm port connecting into the rDNA tangle."""
        ports = self.find_rdna_entry_ports()
        return ports[0] if len(ports) > 0 else None

    def _find_paths(self, entry_port: SegmentPort) -> tuple[list[list[SegmentPort]], list[list[SegmentPort]]]:
        """Search for valid paths from entry_port through rDNA unitigs to distal caps."""
        max_hops = len(self.graph.segments)
        cap_paths: list[list[SegmentPort]] = []
        rdna_paths: list[list[SegmentPort]] = []
        rdna_segments = self._get_rdna_segments()

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
                if self._is_contaminant(nxt.name):
                    continue

                is_cap = self._is_telomeric_cap(nxt) and self._get_usco_count(nxt.name) == 0
                if is_cap:
                    # Distal telomeric cap terminates an rDNA cassette only if rDNA was traversed
                    if any(p.name in rdna_segments for p in path) or (nxt.name in rdna_segments):
                        new_path = list(path)
                        new_path.append(nxt)
                        cap_paths.append(new_path)
                        advanced = True
                    continue

                if nxt.name in rdna_segments:
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

        def score_path(p: list[SegmentPort]) -> tuple[int, int, tuple[str, ...]]:
            rdna_count = sum(1 for port in p if self.is_rdna_unitig(port))
            return (rdna_count, self.compute_path_length(p), tuple(port.to_symbol() for port in p))

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

        if resolved_entry.name in self._get_rdna_segments():
            return [resolved_entry] if include_anchor else []

        cap_paths, rdna_paths = self._find_paths(resolved_entry)

        best_path: list[SegmentPort] = []
        if len(cap_paths) > 0:

            def score_cap(p: list[SegmentPort]) -> tuple[int, int, tuple[str, ...]]:
                cnt = sum(1 for port in p if self.is_rdna_unitig(port))
                return (cnt, self.compute_path_length(p), tuple(port.to_symbol() for port in p))

            cap_paths.sort(key=score_cap, reverse=True)
            chosen = cap_paths[0]
            if include_anchor:
                best_path = chosen
            else:
                best_path = [p for p in chosen if self.is_rdna_unitig(p)]
        elif len(rdna_paths) > 0:

            def score_rdna(p: list[SegmentPort]) -> tuple[int, int, tuple[str, ...]]:
                cnt = sum(1 for port in p if self.is_rdna_unitig(port))
                return (cnt, self.compute_path_length(p), tuple(port.to_symbol() for port in p))

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
                        if edge.cigar and edge.cigar != "*":
                            overlap = cigar_query_consumption(edge.cigar)
                        else:
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

    def resolve_entry_stem(
        self,
        entry_port: SegmentPort,
        chrom: str = "II",
        side: str = "5p",
        gap_index: int = 1,
        notes: str = "",
    ) -> RdnaResolution:
        """Resolve rDNA cassette walk from a core entry stem to a distal cap (graph-connected or heuristic)."""
        path_to_cap = self.reach_distal_telomere_cap(entry_port)

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
            else:
                ev = self._lookup_evidence(path_to_cap[0])
                if ev is not None and ev.depth > 0:
                    core_cov = ev.depth

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
                distal_cap_port=cap_port,
                capped_via_heuristic=False,
                reached_distal_cap=True,
            )

        # Distal cap not reached
        feature_id_end = build_feature_id(EndState.END_SHARED_EXEMPT, chrom, side)
        feature_id_quirk = build_feature_id(IntervalQuirk.RDNA_CASSETTE, chrom, side, index=gap_index)

        unresolved_path = self.traverse_rdna_cassette(entry_port, include_anchor=True)
        all_unitigs = tuple(p.name for p in unresolved_path)
        path_len = self.compute_path_length(unresolved_path) if len(unresolved_path) > 0 else None
        assembled_seq = self.assemble_contiguous_path(unresolved_path) if len(unresolved_path) > 0 else None

        left_str = entry_port.to_symbol()

        distal_cap_port: SegmentPort | None = None
        capped_via_heuristic = False
        chained = [p for p in unresolved_path if self.is_rdna_unitig(p)]
        if len(chained) >= 1:
            expected_dir = chained[0].strand
            visited = {p.name for p in unresolved_path}
            distal_cap_port = self.find_heuristic_rdna_cap(expected_dir, visited)
            if distal_cap_port is not None:
                capped_via_heuristic = True

        if capped_via_heuristic and distal_cap_port is not None:
            right_str = distal_cap_port.to_symbol()
            agp_component = "U"
            agp_evidence = "rna_tangle"
            unresolved_quirk_notes = (
                f"rDNA cassette tangle capped via heuristic at distal cap {distal_cap_port.to_symbol()}"
            )
            unresolved_end_notes = f"Shared exempt rDNA cassette capped via heuristic at {distal_cap_port.to_symbol()}"
        else:
            right_str = None
            agp_component = "U"
            agp_evidence = "unresolved-tangle"
            unresolved_quirk_notes = "Unresolved rDNA cassette tangle: distal telomere cap not reached"
            unresolved_end_notes = "Unresolved rDNA cassette end"

        if notes:
            unresolved_quirk_notes = f"{unresolved_quirk_notes}; {notes}"
            unresolved_end_notes = f"{unresolved_end_notes}; {notes}"

        unresolved_cov: float | None = None
        anchor_seg = self.graph.segments.get(entry_port.name)
        if anchor_seg is not None and anchor_seg.depth > 0:
            unresolved_cov = anchor_seg.depth
        else:
            ev = self._lookup_evidence(entry_port)
            if ev is not None and ev.depth > 0:
                unresolved_cov = ev.depth

        quirk_record = QuirkRecord(
            feature_id=feature_id_quirk,
            chrom=chrom,
            side=side,
            quirk_type=IntervalQuirk.RDNA_CASSETTE,
            length_min_bp=path_len,
            length_median_bp=path_len,
            length_max_bp=path_len,
            left_anchor_port=left_str,
            right_anchor_port=right_str,
            agp_component_type=agp_component,
            agp_linkage_evidence=agp_evidence,
            unitig_ids=all_unitigs,
            core_coverage_bp_median=unresolved_cov,
            representative_sequence_selected=assembled_seq is not None,
            notes=unresolved_quirk_notes,
        )

        end_state_record = EndStateRecord(
            feature_id=feature_id_end,
            chrom=chrom,
            side=side,
            state=EndState.END_SHARED_EXEMPT,
            shared_port=right_str,
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
            distal_cap_port=distal_cap_port,
            capped_via_heuristic=capped_via_heuristic,
            reached_distal_cap=False,
        )

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

        if resolved_entry is None:
            feature_id_end = build_feature_id(EndState.END_SHARED_EXEMPT, chrom, side)
            feature_id_quirk = build_feature_id(IntervalQuirk.RDNA_CASSETTE, chrom, side, index=gap_index)
            quirk_record = QuirkRecord(
                feature_id=feature_id_quirk,
                chrom=chrom,
                side=side,
                quirk_type=IntervalQuirk.RDNA_CASSETTE,
                length_min_bp=None,
                length_median_bp=None,
                length_max_bp=None,
                left_anchor_port=None,
                right_anchor_port=None,
                agp_component_type="U",
                agp_linkage_evidence="unresolved-tangle",
                unitig_ids=(),
                core_coverage_bp_median=None,
                representative_sequence_selected=False,
                notes="Unresolved rDNA cassette tangle: no entry anchor found" if not notes else notes,
            )
            end_state_record = EndStateRecord(
                feature_id=feature_id_end,
                chrom=chrom,
                side=side,
                state=EndState.END_SHARED_EXEMPT,
                shared_port=None,
                resolved=False,
                notes="Unresolved rDNA cassette end: no entry anchor found" if not notes else notes,
            )
            return RdnaResolution(
                chrom=chrom,
                side=side,
                state=EndState.END_SHARED_EXEMPT,
                end_state_record=end_state_record,
                quirk_record=quirk_record,
                path=[],
                sequence=None,
                distal_cap_port=None,
                capped_via_heuristic=False,
                reached_distal_cap=False,
            )

        return self.resolve_entry_stem(
            entry_port=resolved_entry,
            chrom=chrom,
            side=side,
            gap_index=gap_index,
            notes=notes,
        )

    def collect_unresolved_tangles(self, chained_segments: set[str] | Sequence[str] = ()) -> list[QuirkRecord]:
        """Record any remaining eukaryotic rDNA segments not reached from an entry stem as single-segment unresolved-tangle quirks."""
        chained = set(chained_segments)
        rdna_segs = sorted(self._get_rdna_segments())
        quirks: list[QuirkRecord] = []
        for seg_name in rdna_segs:
            if seg_name in chained:
                continue
            seg = self.graph.segments.get(seg_name)
            seg_len = seg.length if seg is not None else None
            q_fid = build_feature_id(IntervalQuirk.RDNA_CASSETTE, seg_name, "mid", index=1)
            q_rec = QuirkRecord(
                feature_id=q_fid,
                chrom=seg_name,
                side="mid",
                quirk_type=IntervalQuirk.RDNA_CASSETTE,
                length_min_bp=seg_len,
                length_median_bp=seg_len,
                length_max_bp=seg_len,
                left_anchor_port=f"{seg_name}+",
                right_anchor_port=None,
                agp_component_type="U",
                agp_linkage_evidence="unresolved-tangle",
                unitig_ids=(seg_name,),
                core_coverage_bp_median=seg.depth if seg is not None else None,
                representative_sequence_selected=False,
                notes=f"Unresolved unanchored eukaryotic rDNA segment: {seg_name}",
            )
            quirks.append(q_rec)
        return quirks
