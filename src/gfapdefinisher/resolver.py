"""ONT spanning-read resolution engine for variable satellite gaps and orphan caps.

Implements BIO-001 (unconditional terminal completion) and BIO-004 (representative
repeat preservation for variable-length tandem satellites). Bridges candidate unitig
ports using long ONT reads, detects innermost boundary hits to disambiguate multivalent
satellite repeat arrays, computes empirical gap distributions, and emits ratified
QuirkRecord and EndStateRecord instances.
"""

from __future__ import annotations

import gzip
import statistics
from collections import Counter, defaultdict
from collections.abc import Iterable, Iterator, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, TextIO

from gfapdefinisher.emitter import DEFAULT_AGP_POLICY, AgpPolicy
from gfapdefinisher.models import SegmentPort, Strand
from gfapdefinisher.policy import DEFAULT_POLICY, Policy
from gfapdefinisher.quirks import (
    EndState,
    EndStateRecord,
    IntervalQuirk,
    QuirkRecord,
    build_feature_id,
)

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of a nucleotide sequence."""
    return sequence.translate(_COMPLEMENT)[::-1]


def _parse_port_spec(port: str | SegmentPort) -> tuple[str, Strand]:
    """Extract contig name and port strand from a port specification.

    Handles SegmentPort instances, string representations like 'utg000056l+' or 'utg000126l-',
    and plain contig names (defaulting to Strand.FORWARD).
    """
    if isinstance(port, SegmentPort):
        return port.name, port.strand
    s = str(port).strip()
    if s.endswith("+"):
        return s[:-1], Strand.FORWARD
    if s.endswith("-"):
        return s[:-1], Strand.REVERSE
    return s, Strand.FORWARD


def _port_boundary_slack(hit: PafHit, is_left: bool, strand: Strand) -> int:
    """Compute unaligned distance between alignment and segment port boundary.

    For left port (walk exiting port at 3'):
      - '+' strand: boundary is target tlen -> slack is tlen - tend (left_slack).
      - '-' strand: boundary is target 0 -> slack is tstart (right_slack).

    For right port (walk entering port at 5'):
      - '+' strand: boundary is target 0 -> slack is tstart (right_slack).
      - '-' strand: boundary is target tlen -> slack is tlen - tend (left_slack).
    """
    if is_left:
        return hit.left_slack if strand == Strand.FORWARD else hit.right_slack
    return hit.right_slack if strand == Strand.FORWARD else hit.left_slack


@dataclass(slots=True, frozen=True)
class PafHit:
    """A single pairwise mapping format (PAF) alignment hit.

    Attributes:
        qname: Query read identifier.
        qlen: Total query sequence length.
        qstart: Query 0-based alignment start.
        qend: Query 0-based alignment end.
        strand: Alignment relative strand ('+' or '-').
        tname: Target sequence identifier.
        tlen: Total target sequence length.
        tstart: Target 0-based alignment start.
        tend: Target 0-based alignment end.
        nmatch: Number of matching residue bases.
        alen: Alignment block length.
        mapq: Mapping quality (0-255).
        tags: SAM-like key-value optional tags.
    """

    qname: str
    qlen: int
    qstart: int
    qend: int
    strand: str
    tname: str
    tlen: int
    tstart: int
    tend: int
    nmatch: int
    alen: int
    mapq: int
    tags: dict[str, str] = field(default_factory=dict)

    @classmethod
    def from_line(cls, line: str) -> PafHit:
        """Parse a single PAF format line into a PafHit instance.

        Format column indices (0..11) adhere to standard PAF format specifications.
        """
        parts = line.rstrip("\r\n").split("\t")
        if len(parts) < 12:
            raise ValueError(f"Malformed PAF line: expected >= 12 columns, got {len(parts)}")

        tags: dict[str, str] = {}
        for tag_str in parts[12:]:
            if ":" in tag_str:
                tag_name = tag_str.split(":", 1)[0]
                tags[tag_name] = tag_str

        return cls(
            qname=parts[0],
            qlen=int(parts[1]),
            qstart=int(parts[2]),
            qend=int(parts[3]),
            strand=parts[4],
            tname=parts[5],
            tlen=int(parts[6]),
            tstart=int(parts[7]),
            tend=int(parts[8]),
            nmatch=int(parts[9]),
            alen=int(parts[10]),
            mapq=int(parts[11]),
            tags=tags,
        )

    @property
    def left_slack(self) -> int:
        """Unaligned distance from alignment end to target 3' terminus."""
        return max(0, self.tlen - self.tend)

    @property
    def right_slack(self) -> int:
        """Unaligned distance from target 5' terminus to alignment start."""
        return max(0, self.tstart)


def parse_paf_hits(paf_file: Path | str | TextIO | Iterable[str]) -> Iterator[PafHit]:
    """Stream PAF records lazily with minimal memory footprint.

    Args:
        paf_file: Input file path, open text stream, or line iterable.

    Yields:
        Parsed PafHit records line by line.
    """
    if isinstance(paf_file, str | Path):
        path = Path(paf_file)
        if str(path).endswith(".gz"):
            with gzip.open(path, "rt", encoding="utf-8") as handle:
                for line in handle:
                    s = line.strip()
                    if s and not s.startswith("#"):
                        yield PafHit.from_line(s)
        else:
            with open(path, encoding="utf-8") as handle:
                for line in handle:
                    s = line.strip()
                    if s and not s.startswith("#"):
                        yield PafHit.from_line(s)
    else:
        for line in paf_file:
            s = line.strip()
            if s and not s.startswith("#"):
                yield PafHit.from_line(s)


@dataclass(slots=True, frozen=True)
class SpanningRead:
    """A long sequencing read bridging between candidate left and right flank ports.

    Attributes:
        qname: Query read identifier.
        strand: Orientation strand relative to flank walk ('+' or '-').
        left_hit: Innermost alignment hit to left flank port.
        right_hit: Innermost alignment hit to right flank port.
        gap_length: Spliced gap interval length in base pairs.
        q_from: Query start coordinate of gap interval.
        q_to: Query end coordinate of gap interval.
        sequence: Extracted representative gap sequence if baited or provided.
    """

    qname: str
    strand: str
    left_hit: PafHit
    right_hit: PafHit
    gap_length: int
    q_from: int
    q_to: int
    sequence: str | None = None


@dataclass(slots=True, frozen=True)
class GapDistribution:
    """Empirical distribution metrics across long spanning reads."""

    min_bp: int
    median_bp: int
    max_bp: int
    iqr_bp: int
    count: int


@dataclass(slots=True)
class SatelliteGapResolution:
    """Typed resolution artifact for a variable tandem satellite gap locus."""

    chrom: str
    side: str
    state: EndState
    quirk_record: QuirkRecord
    end_state_record: EndStateRecord
    spanning_reads: list[SpanningRead]
    representative_read: SpanningRead | None = None
    representative_sequence: str | None = None
    distribution: GapDistribution | None = None

    @property
    def left_flank_port(self) -> SegmentPort | None:
        if self.quirk_record and self.quirk_record.left_anchor_port:
            name, strand = _parse_port_spec(self.quirk_record.left_anchor_port)
            return SegmentPort(name, strand)
        return None

    @property
    def right_flank_port(self) -> SegmentPort | None:
        if self.quirk_record and self.quirk_record.right_anchor_port:
            name, strand = _parse_port_spec(self.quirk_record.right_anchor_port)
            return SegmentPort(name, strand)
        return None

    @property
    def patched_sequence(self) -> str | None:
        return self.representative_sequence

    @property
    def median_gap_bp(self) -> float:
        if self.quirk_record and self.quirk_record.length_median_bp is not None:
            return float(self.quirk_record.length_median_bp)
        return float(DEFAULT_AGP_POLICY.unknown_gap_bp.value)

    @property
    def gap_index(self) -> int:
        if self.quirk_record and self.quirk_record.feature_id:
            parts = self.quirk_record.feature_id.split(":")
            if len(parts) > 2 and parts[-1].isdigit():
                return int(parts[-1])
        return 1


GapResolution = SatelliteGapResolution


def load_reads_sequences(
    reads_path: Path | str | None,
    target_qnames: set[str] | Sequence[str] | None = None,
) -> dict[str, str]:
    """Load long reads into memory as a mapping from read name to sequence string.

    Supports FASTA and FASTQ formats, uncompressed or gzip-compressed.
    If target_qnames is provided, only reads whose IDs are in target_qnames are retained,
    enabling memory-efficient single-pass streaming across large FASTQ/FASTA files.
    """
    if reads_path is None:
        return {}

    path = Path(reads_path)
    if not path.is_file():
        return {}

    target_set = set(target_qnames) if target_qnames is not None else None
    if target_set is not None and len(target_set) == 0:
        return {}

    is_gz = path.suffix in (".gz", ".gzip")
    open_fn = gzip.open if is_gz else open

    sequences: dict[str, str] = {}
    with open_fn(path, "rt", encoding="utf-8") as f:
        first_char = ""
        for line in f:
            stripped = line.strip()
            if stripped:
                first_char = stripped[0]
                break
        f.seek(0)

        if first_char == ">":
            curr_id = ""
            curr_seq: list[str] = []
            capture = False
            for line in f:
                line_str = line.strip()
                if not line_str:
                    continue
                if line_str.startswith(">"):
                    if curr_id and capture:
                        sequences[curr_id] = "".join(curr_seq)
                        if target_set is not None and len(sequences) >= len(target_set):
                            break
                    curr_id = line_str[1:].split()[0]
                    capture = target_set is None or curr_id in target_set
                    curr_seq = []
                elif capture:
                    curr_seq.append(line_str)
            if curr_id and capture and (target_set is None or curr_id in target_set):
                sequences[curr_id] = "".join(curr_seq)
        elif first_char == "@":
            while True:
                header = f.readline()
                if not header:
                    break
                header_str = header.strip()
                if not header_str:
                    continue
                qname = header_str[1:].split()[0]
                seq = f.readline().strip()
                f.readline()  # + separator
                f.readline()  # qual
                if target_set is None or qname in target_set:
                    sequences[qname] = seq
                    if target_set is not None and len(sequences) >= len(target_set):
                        break

    return sequences


class OntBridgingEngine:
    """Resolution engine bridging assembly graph ports with long ONT reads.

    Provides innermost boundary hit detection to resolve multivalent tandem satellite
    repeat alignments, empirical gap length distribution statistics, and BIO-004
    representative sequence selection.
    """

    def __init__(
        self,
        graph_or_policy: Any = DEFAULT_POLICY,
        ledger: Any | None = None,
        policy: Policy = DEFAULT_POLICY,
        *,
        graph: Any | None = None,
    ) -> None:
        if graph is not None:
            self.graph = graph
            self.ledger = ledger
            self.policy = graph_or_policy if isinstance(graph_or_policy, Policy) else policy
        elif isinstance(graph_or_policy, Policy):
            self.policy = graph_or_policy
            self.graph = None
            self.ledger = ledger
        elif hasattr(graph_or_policy, "segments"):
            self.graph = graph_or_policy
            self.ledger = ledger
            self.policy = policy
        else:
            self.graph = None
            self.ledger = ledger
            self.policy = policy

        self.gap_resolutions: list[SatelliteGapResolution] = []
        self.dropout_resolutions: list[tuple[QuirkRecord, EndStateRecord, str | None]] = []

    def collect_spanning_read_ids(
        self,
        paf_source: Path | str | TextIO | Iterable[str] | Iterable[PafHit],
    ) -> set[str]:
        """Collect read IDs (qnames) that bridge candidate port pairs."""
        source_iterable: Iterable[str | PafHit]
        if isinstance(paf_source, str | Path | TextIO):
            source_iterable = parse_paf_hits(paf_source)
        else:
            source_iterable = paf_source

        all_hits: list[PafHit] = []
        hits_by_read: dict[str, list[PafHit]] = defaultdict(list)
        for item in source_iterable:
            if isinstance(item, PafHit):
                hit = item
            else:
                s = item.strip()
                if not s or s.startswith("#"):
                    continue
                hit = PafHit.from_line(s)
            all_hits.append(hit)
            hits_by_read[hit.qname].append(hit)

        pair_counts: dict[tuple[str, str], int] = defaultdict(int)
        for read_hits in hits_by_read.values():
            targets = sorted({h.tname for h in read_hits})
            for i in range(len(targets)):
                for j in range(i + 1, len(targets)):
                    pair_counts[(targets[i], targets[j])] += 1

        target_qnames: set[str] = set()
        checked_pairs: set[tuple[str, str]] = set()
        for t1, t2 in pair_counts:
            if (t1, t2) in checked_pairs or (t2, t1) in checked_pairs:
                continue
            checked_pairs.add((t1, t2))
            spanning = self.collect_spanning_reads(paf_source=all_hits, left_port=t1, right_port=t2)
            if not spanning:
                spanning = self.collect_spanning_reads(paf_source=all_hits, left_port=t2, right_port=t1)
            if len(spanning) >= self.policy.ont.min_spanning_reads.value:
                for sr in spanning:
                    target_qnames.add(sr.qname)

        return target_qnames

    def identify_orphan_caps(self) -> dict[str, tuple[str, str]]:
        """Identify orphan caps in graph.

        Returns a mapping from segment_name to (telomere_port, inward_port).
        """
        orphan_caps: dict[str, tuple[str, str]] = {}
        if self.graph is None or not self.ledger:
            return orphan_caps

        ledger_map = self.ledger
        if isinstance(ledger_map, list | tuple):
            ledger_map = {ev.port.to_symbol(): ev for ev in ledger_map if hasattr(ev, "port")}

        for seg_name in self.graph.segments:
            deg_plus = self.graph.degree(SegmentPort(seg_name, Strand.FORWARD))
            deg_minus = self.graph.degree(SegmentPort(seg_name, Strand.MINUS))
            if deg_plus != 0 or deg_minus != 0:
                continue

            ev_plus = ledger_map.get(f"{seg_name}+")
            ev_minus = ledger_map.get(f"{seg_name}-")

            if getattr(ev_plus, "usco_count", 0) != 0 or getattr(ev_minus, "usco_count", 0) != 0:
                continue
            if getattr(ev_plus, "usco_hit_count", 0) != 0 or getattr(ev_minus, "usco_hit_count", 0) != 0:
                continue

            is_contam_plus = getattr(ev_plus, "is_contaminant", False)
            if callable(is_contam_plus):
                is_contam_plus = is_contam_plus()
            is_contam_minus = getattr(ev_minus, "is_contaminant", False)
            if callable(is_contam_minus):
                is_contam_minus = is_contam_minus()
            if is_contam_plus or is_contam_minus:
                continue

            has_plus = ev_plus is not None and getattr(ev_plus, "has_terminal_telomere", False)
            has_minus = ev_minus is not None and getattr(ev_minus, "has_terminal_telomere", False)
            if has_minus == has_plus:
                continue

            if has_minus:
                orphan_caps[seg_name] = (f"{seg_name}-", f"{seg_name}+")
            elif has_plus:
                orphan_caps[seg_name] = (f"{seg_name}+", f"{seg_name}-")

        return orphan_caps

    def identify_unlinked_grs_unitigs(self) -> set[str]:
        """Identify unlinked GRS repeat unitigs in graph.

        Returns segment IDs that have degree 0 at both ports, zero USCOs,
        non-contaminant status, no rDNA hits, and no terminal telomeres.
        """
        unlinked_grs: set[str] = set()
        if self.graph is None or not self.ledger:
            return unlinked_grs

        ledger_map = self.ledger
        if isinstance(ledger_map, list | tuple):
            ledger_map = {ev.port.to_symbol(): ev for ev in ledger_map if hasattr(ev, "port")}

        for seg_name in self.graph.segments:
            deg_plus = self.graph.degree(SegmentPort(seg_name, Strand.FORWARD))
            deg_minus = self.graph.degree(SegmentPort(seg_name, Strand.MINUS))
            if deg_plus != 0 or deg_minus != 0:
                continue

            ev_plus = ledger_map.get(f"{seg_name}+")
            ev_minus = ledger_map.get(f"{seg_name}-")

            if getattr(ev_plus, "usco_count", 0) != 0 or getattr(ev_minus, "usco_count", 0) != 0:
                continue
            if getattr(ev_plus, "usco_hit_count", 0) != 0 or getattr(ev_minus, "usco_hit_count", 0) != 0:
                continue

            if getattr(ev_plus, "rdna_hit_count", 0) != 0 or getattr(ev_minus, "rdna_hit_count", 0) != 0:
                continue
            if getattr(ev_plus, "has_terminal_rdna", False) or getattr(ev_minus, "has_terminal_rdna", False):
                continue

            is_contam_plus = getattr(ev_plus, "is_contaminant", False)
            if callable(is_contam_plus):
                is_contam_plus = is_contam_plus()
            is_contam_minus = getattr(ev_minus, "is_contaminant", False)
            if callable(is_contam_minus):
                is_contam_minus = is_contam_minus()
            if is_contam_plus or is_contam_minus:
                continue

            has_plus = ev_plus is not None and getattr(ev_plus, "has_terminal_telomere", False)
            has_minus = ev_minus is not None and getattr(ev_minus, "has_terminal_telomere", False)
            if has_plus or has_minus:
                continue

            unlinked_grs.add(seg_name)

        return unlinked_grs

    def _has_terminal_telomere(self, port: SegmentPort | str) -> bool:
        """Check whether a port carries a terminal telomere in the evidence ledger."""
        if not self.ledger:
            return False
        port_name, port_strand = _parse_port_spec(port)
        port_sym = f"{port_name}{port_strand.to_symbol()}"
        ledger_map = self.ledger
        if isinstance(ledger_map, list | tuple):
            ledger_map = {ev.port.to_symbol(): ev for ev in ledger_map if hasattr(ev, "port")}
        ev = ledger_map.get(port_sym)
        return bool(ev is not None and getattr(ev, "has_terminal_telomere", False))

    def identify_and_resolve_gaps(
        self,
        paf_source: Path | str | TextIO | Iterable[str] | Iterable[PafHit],
        end_states: list[EndStateRecord],
        quirks: list[QuirkRecord],
        reads_sequences: Mapping[str, str] | Path | str | None = None,
    ) -> tuple[list[EndStateRecord], list[QuirkRecord]]:
        """Identify candidate variable satellite gap pairs and resolve them via ONT bridging."""
        orphan_caps = self.identify_orphan_caps()
        unlinked_grs = self.identify_unlinked_grs_unitigs()
        resolved_end_keys: set[tuple[str, str]] = set()
        for er in end_states:
            if er.resolved:
                resolved_end_keys.add((er.chrom, er.side))
                if hasattr(er, "segment_id"):
                    resolved_end_keys.add((er.segment_id, er.side))

        hits_by_read: dict[str, list[PafHit]] = defaultdict(list)
        source_iterable: Iterable[str | PafHit]
        if isinstance(paf_source, str | Path | TextIO):
            source_iterable = parse_paf_hits(paf_source)
        else:
            source_iterable = paf_source

        all_hits: list[PafHit] = []
        for item in source_iterable:
            if isinstance(item, PafHit):
                hit = item
            else:
                s = item.strip()
                if not s or s.startswith("#"):
                    continue
                hit = PafHit.from_line(s)
            all_hits.append(hit)
            hits_by_read[hit.qname].append(hit)

        active_reads: Mapping[str, str] | None = None
        if isinstance(reads_sequences, str | Path):
            target_ids = self.collect_spanning_read_ids(all_hits)
            active_reads = load_reads_sequences(reads_sequences, target_qnames=target_ids)
        else:
            active_reads = reads_sequences

        pair_counts: dict[tuple[str, str], int] = defaultdict(int)
        for read_hits in hits_by_read.values():
            targets = sorted({h.tname for h in read_hits})
            for i in range(len(targets)):
                for j in range(i + 1, len(targets)):
                    pair_counts[(targets[i], targets[j])] += 1

        sorted_pairs = sorted(
            pair_counts.items(),
            key=lambda kv: (-kv[1], kv[0][0], kv[0][1]),
        )

        resolved_pairs: set[tuple[str, str]] = set()
        used_cap_segs: set[str] = set()

        for (t1, t2), _count in sorted_pairs:
            if (t1, t2) in resolved_pairs or (t2, t1) in resolved_pairs:
                continue

            matched_end: EndStateRecord | None = None
            partner_target: str = ""

            for er in end_states:
                if er.resolved or er.germline_port is None:
                    continue
                end_key = (getattr(er, "segment_id", er.chrom), er.side)
                if end_key in resolved_end_keys or (er.chrom, er.side) in resolved_end_keys:
                    continue

                anchor = er.germline_port
                anchor_clean = anchor.rstrip("+-")
                t1_clean = t1.rstrip("+-")
                t2_clean = t2.rstrip("+-")

                t1_matches = t1 == anchor or t1_clean == anchor_clean
                t2_matches = t2 == anchor or t2_clean == anchor_clean

                if not t1_matches and not t2_matches:
                    continue

                if t1_matches and not t2_matches:
                    partner_target = t2
                elif t2_matches and not t1_matches:
                    partner_target = t1
                else:
                    partner_target = t2

                partner_seg = partner_target.rstrip("+-")
                if orphan_caps:
                    if partner_seg not in orphan_caps or partner_seg in used_cap_segs:
                        continue

                matched_end = er
                break

            if matched_end is None:
                if end_states:
                    continue
                cand_left = t1
                cand_right = t2

            best_spanning: list[SpanningRead] = []
            best_left: SegmentPort | str = t1
            best_right: SegmentPort | str = t2

            if matched_end is not None and orphan_caps and partner_target:
                partner_seg = partner_target.rstrip("+-")
                anchor_name, anchor_strand = _parse_port_spec(matched_end.germline_port)
                anchor_ports = [SegmentPort(anchor_name, anchor_strand)]

                cap_tel_port, cap_inw_port = orphan_caps[partner_seg]
                cap_name, tel_strand = _parse_port_spec(cap_tel_port)
                cap_ports = [SegmentPort(cap_name, tel_strand), SegmentPort(cap_name, tel_strand.invert())]

                for a_port in anchor_ports:
                    for c_port in cap_ports:
                        sp = self.collect_spanning_reads(
                            paf_source=all_hits,
                            left_port=a_port,
                            right_port=c_port,
                            reads_sequences=active_reads,
                        )
                        if len(sp) > len(best_spanning):
                            best_spanning = sp
                            best_left = a_port
                            best_right = c_port

                spanning = best_spanning
                left_port, right_port = best_left, best_right
            else:
                left_port, right_port = cand_left, cand_right
                spanning = self.collect_spanning_reads(
                    paf_source=all_hits,
                    left_port=cand_left,
                    right_port=cand_right,
                    reads_sequences=active_reads,
                )
                if not spanning:
                    spanning = self.collect_spanning_reads(
                        paf_source=all_hits,
                        left_port=cand_right,
                        right_port=cand_left,
                        reads_sequences=active_reads,
                    )
                    if spanning and (matched_end is None or not orphan_caps):
                        left_port, right_port = cand_right, cand_left

            floor = self.policy.ont.min_spanning_reads.value
            if len(spanning) >= floor:
                resolved_pairs.add((t1, t2))
                resolved_pairs.add((t2, t1))
                resolved_pairs.add((str(left_port), str(right_port)))
                resolved_pairs.add((str(right_port), str(left_port)))
                if matched_end is not None and partner_target:
                    used_cap_segs.add(partner_target.rstrip("+-"))

                chrom = matched_end.chrom if matched_end else str(left_port).rstrip("+-")
                side = matched_end.side if matched_end else ("3p" if "+" in str(left_port) else "5p")

                resolved_end_keys.add((chrom, side))
                if matched_end and hasattr(matched_end, "segment_id"):
                    resolved_end_keys.add((matched_end.segment_id, side))

                res = self.resolve_satellite_gap(
                    chrom=chrom,
                    side=side,
                    left_port=left_port,
                    right_port=right_port,
                    paf_source=all_hits,
                    spanning_reads=spanning,
                    reads_sequences=active_reads,
                    core_coverage_bp_median=matched_end.core_coverage_bp_median if matched_end else None,
                    core_depth_ratio=matched_end.germline_depth_ratio if matched_end else None,
                    somatic_depth_ratio=matched_end.somatic_depth_ratio if matched_end else None,
                    somatic_port=matched_end.somatic_port if matched_end else None,
                    matched_end=matched_end,
                )

                new_ends: list[EndStateRecord] = []
                replaced_end = False
                for er in end_states:
                    if er.chrom == chrom and er.side == side:
                        new_ends.append(res.end_state_record)
                        replaced_end = True
                    else:
                        new_ends.append(er)
                if not replaced_end:
                    new_ends.append(res.end_state_record)
                end_states = new_ends

                new_quirks: list[QuirkRecord] = []
                replaced_quirk = False
                for qr in quirks:
                    if qr.chrom == chrom and qr.side == side and qr.quirk_type == IntervalQuirk.SAT_GAP_VARIABLE:
                        new_quirks.append(res.quirk_record)
                        replaced_quirk = True
                    else:
                        new_quirks.append(qr)
                if not replaced_quirk:
                    new_quirks.append(res.quirk_record)
                quirks = new_quirks

        floor = self.policy.ont.min_spanning_reads.value
        for er in end_states:
            end_key = (getattr(er, "segment_id", er.chrom), er.side)
            if (
                er.resolved
                or er.germline_port is None
                or end_key in resolved_end_keys
                or (er.chrom, er.side) in resolved_end_keys
            ):
                continue

            anchor_name, anchor_strand = _parse_port_spec(er.germline_port)
            anchor_port = SegmentPort(anchor_name, anchor_strand)

            candidate_chains: list[tuple[int, str, str]] = []
            for mid_id in unlinked_grs:
                cnt1 = pair_counts.get((anchor_name, mid_id), 0) + pair_counts.get((mid_id, anchor_name), 0)
                if cnt1 < floor:
                    continue
                for cap_id in orphan_caps:
                    if cap_id in used_cap_segs:
                        continue
                    cnt2 = pair_counts.get((mid_id, cap_id), 0) + pair_counts.get((cap_id, mid_id), 0)
                    if cnt2 < floor:
                        continue
                    candidate_chains.append((cnt1 + cnt2, mid_id, cap_id))

            candidate_chains.sort(key=lambda item: (-item[0], item[1], item[2]))

            for _score, mid_id, cap_id in candidate_chains:
                if cap_id in used_cap_segs:
                    continue

                best_hop1: tuple[SegmentPort, list[SpanningRead]] | None = None
                for mid_strand in (Strand.FORWARD, Strand.REVERSE):
                    mid_entry = SegmentPort(mid_id, mid_strand)
                    sp1 = self.collect_spanning_reads(
                        paf_source=all_hits,
                        left_port=anchor_port,
                        right_port=mid_entry,
                        reads_sequences=active_reads,
                    )
                    if len(sp1) >= floor:
                        if best_hop1 is None or len(sp1) > len(best_hop1[1]):
                            best_hop1 = (mid_entry, sp1)

                if best_hop1 is None:
                    continue

                mid_entry, sp1 = best_hop1
                mid_exit = mid_entry.invert()

                cap_tel_port, cap_inw_port = orphan_caps[cap_id]
                cap_name, tel_strand = _parse_port_spec(cap_tel_port)
                cap_ports = [SegmentPort(cap_name, tel_strand), SegmentPort(cap_name, tel_strand.invert())]

                best_hop2: tuple[SegmentPort, list[SpanningRead]] | None = None
                for c_port in cap_ports:
                    sp2 = self.collect_spanning_reads(
                        paf_source=all_hits,
                        left_port=mid_exit,
                        right_port=c_port,
                        reads_sequences=active_reads,
                    )
                    if len(sp2) >= floor and self._has_terminal_telomere(c_port):
                        if best_hop2 is None or len(sp2) > len(best_hop2[1]):
                            best_hop2 = (c_port, sp2)

                if best_hop2 is None:
                    continue

                c_port, sp2 = best_hop2

                res1 = self.resolve_satellite_gap(
                    chrom=er.chrom,
                    side=er.side,
                    left_port=anchor_port,
                    right_port=mid_entry,
                    paf_source=all_hits,
                    gap_index=1,
                    spanning_reads=sp1,
                    reads_sequences=active_reads,
                    core_coverage_bp_median=er.core_coverage_bp_median,
                    core_depth_ratio=er.germline_depth_ratio,
                    somatic_depth_ratio=er.somatic_depth_ratio,
                    somatic_port=er.somatic_port,
                    matched_end=er,
                )

                res2 = self.resolve_satellite_gap(
                    chrom=er.chrom,
                    side=er.side,
                    left_port=mid_exit,
                    right_port=c_port,
                    paf_source=all_hits,
                    gap_index=2,
                    spanning_reads=sp2,
                    reads_sequences=active_reads,
                    core_coverage_bp_median=er.core_coverage_bp_median,
                    core_depth_ratio=er.germline_depth_ratio,
                    somatic_depth_ratio=er.somatic_depth_ratio,
                    somatic_port=er.somatic_port,
                    matched_end=er,
                )

                resolved_end_keys.add((er.chrom, er.side))
                if hasattr(er, "segment_id"):
                    resolved_end_keys.add((er.segment_id, er.side))
                used_cap_segs.add(cap_id)

                new_ends = []
                for e_rec in end_states:
                    if e_rec.chrom == er.chrom and e_rec.side == er.side:
                        new_ends.append(res2.end_state_record)
                    else:
                        new_ends.append(e_rec)
                end_states = new_ends

                quirks = [
                    q
                    for q in quirks
                    if not (
                        q.chrom == er.chrom and q.side == er.side and q.quirk_type == IntervalQuirk.SAT_GAP_VARIABLE
                    )
                ]
                quirks.append(res1.quirk_record)
                quirks.append(res2.quirk_record)
                break

        end_states, quirks = self.resolve_subtelomeric_satellite_caps(
            end_states=end_states,
            quirks=quirks,
            orphan_caps=orphan_caps,
            used_cap_segs=used_cap_segs,
            resolved_end_keys=resolved_end_keys,
        )

        return end_states, quirks

    def resolve_subtelomeric_satellite_caps(
        self,
        end_states: list[EndStateRecord],
        quirks: list[QuirkRecord],
        orphan_caps: dict[str, tuple[str, str]] | None = None,
        used_cap_segs: set[str] | None = None,
        resolved_end_keys: set[tuple[str, str]] | None = None,
    ) -> tuple[list[EndStateRecord], list[QuirkRecord]]:
        """Resolve remaining unresolved ends whose distal tip terminates in a tandem satellite repeat matching an orphan cap."""
        if self.graph is None:
            return end_states, quirks
        if orphan_caps is None:
            orphan_caps = self.identify_orphan_caps()
        if used_cap_segs is None:
            used_cap_segs = set()
        if resolved_end_keys is None:
            resolved_end_keys = set()
            for er in end_states:
                if er.resolved:
                    resolved_end_keys.add((er.chrom, er.side))
                    if hasattr(er, "segment_id"):
                        resolved_end_keys.add((er.segment_id, er.side))

        window = self.policy.ont.min_anchor_bp.value * 2
        k = (2 + 2 + 1) * (2 + 1)
        min_repeats = 2 + 2 + 1
        floor = self.policy.ont.min_spanning_reads.value
        unknown_gap = getattr(self.policy, "agp", AgpPolicy()).unknown_gap_bp.value
        canonical_fwd = self.policy.telomere.motif_forward.value
        canonical_rev = self.policy.telomere.motif_reverse.value

        for er in list(end_states):
            end_key = (getattr(er, "segment_id", er.chrom), er.side)
            if (
                er.resolved
                or er.germline_port is None
                or end_key in resolved_end_keys
                or (er.chrom, er.side) in resolved_end_keys
            ):
                continue

            anchor_name, anchor_strand = _parse_port_spec(er.germline_port)
            anchor_port = SegmentPort(anchor_name, anchor_strand)
            anchor_seg = self.graph.segments.get(anchor_name)
            if anchor_seg is None or not anchor_seg.sequence:
                continue

            anchor_seq = (
                anchor_seg.sequence if anchor_strand == Strand.FORWARD else reverse_complement(anchor_seg.sequence)
            )
            if len(anchor_seq) < window:
                continue

            tip_seq = anchor_seq[-window:]
            kmers = Counter(tip_seq[i : i + k] for i in range(len(tip_seq) - k + 1))
            if not kmers:
                continue
            top_kmer, top_cnt = kmers.most_common(1)[0]
            if top_cnt < min_repeats:
                continue
            if canonical_fwd in top_kmer or canonical_rev in top_kmer:
                continue

            rc_top = reverse_complement(top_kmer)

            for cap_id in sorted(orphan_caps.keys() - used_cap_segs):
                cap_seg = self.graph.segments.get(cap_id)
                if cap_seg is None or not cap_seg.sequence:
                    continue

                cap_seq = cap_seg.sequence
                cnt_in_cap = cap_seq.count(top_kmer) + cap_seq.count(rc_top)
                if cnt_in_cap < floor:
                    continue

                cap_tel_port, _ = orphan_caps[cap_id]
                cap_tel_name, cap_tel_strand = _parse_port_spec(cap_tel_port)
                cap_port = SegmentPort(cap_tel_name, cap_tel_strand)

                feature_id_quirk = build_feature_id(IntervalQuirk.SAT_GAP_VARIABLE, er.chrom, er.side, index=1)
                feature_id_end = build_feature_id(EndState.END_COMPLETE, er.chrom, er.side)

                quirk_record = QuirkRecord(
                    feature_id=feature_id_quirk,
                    chrom=er.chrom,
                    side=er.side,
                    quirk_type=IntervalQuirk.SAT_GAP_VARIABLE,
                    spanning_read_count=0,
                    length_min_bp=unknown_gap,
                    length_median_bp=unknown_gap,
                    length_max_bp=unknown_gap,
                    length_iqr_bp=0,
                    repeat_monomer_len_bp=None,
                    copy_number_range=None,
                    left_anchor_port=anchor_port.to_symbol(),
                    right_anchor_port=cap_port.to_symbol(),
                    agp_component_type="U",
                    agp_linkage_evidence="paired-ends",
                    unitig_ids=(anchor_name, cap_id),
                    read_count=0,
                    core_coverage_bp_median=er.core_coverage_bp_median,
                    core_depth_ratio=er.germline_depth_ratio,
                    representative_sequence_selected=False,
                    sequence_provenance="HiFi-germline",
                    unpolished=False,
                    notes="Resolved via subtelomeric satellite capping heuristic; unbridged repeat gap (type U)",
                )

                end_state_record = EndStateRecord(
                    feature_id=feature_id_end,
                    chrom=er.chrom,
                    side=er.side,
                    state=EndState.END_COMPLETE,
                    somatic_port=er.somatic_port,
                    germline_port=cap_port.to_symbol(),
                    shared_port=None,
                    core_coverage_bp_median=er.core_coverage_bp_median,
                    somatic_depth_ratio=er.somatic_depth_ratio,
                    germline_depth_ratio=er.germline_depth_ratio,
                    resolved=True,
                    notes=f"Resolved via subtelomeric satellite capping heuristic to orphan cap {cap_id}",
                )

                res = SatelliteGapResolution(
                    chrom=er.chrom,
                    side=er.side,
                    state=EndState.END_COMPLETE,
                    quirk_record=quirk_record,
                    end_state_record=end_state_record,
                    spanning_reads=[],
                    representative_read=None,
                    representative_sequence=None,
                    distribution=None,
                )
                self.gap_resolutions.append(res)
                resolved_end_keys.add((er.chrom, er.side))
                if hasattr(er, "segment_id"):
                    resolved_end_keys.add((er.segment_id, er.side))
                used_cap_segs.add(cap_id)

                new_ends = []
                for e_rec in end_states:
                    if e_rec.chrom == er.chrom and e_rec.side == er.side:
                        new_ends.append(end_state_record)
                    else:
                        new_ends.append(e_rec)
                end_states = new_ends

                quirks = [
                    q
                    for q in quirks
                    if not (
                        q.chrom == er.chrom and q.side == er.side and q.quirk_type == IntervalQuirk.SAT_GAP_VARIABLE
                    )
                ]
                quirks.append(quirk_record)
                break

        return end_states, quirks

    def collect_spanning_reads(
        self,
        paf_source: Path | str | TextIO | Iterable[str] | Iterable[PafHit],
        left_port: str | SegmentPort,
        right_port: str | SegmentPort,
        *,
        min_anchor_bp: int | None = None,
        max_boundary_slack: int | None = None,
        min_mapq: int | None = None,
        reads_sequences: Mapping[str, str] | None = None,
    ) -> list[SpanningRead]:
        """Identify long reads bridging candidate ports with innermost boundary hit detection.

        Filters alignments by anchor length and boundary slack, enforces dovetail
        collinearity, and resolves multivalent tandem satellite repeats by selecting
        the innermost boundary hits closest to the intervening gap.

        Args:
            paf_source: PAF file path, text stream, line iterable, or PafHit iterable.
            left_port: Left flank target identifier or SegmentPort.
            right_port: Right flank target identifier or SegmentPort.
            min_anchor_bp: Minimum alignment length threshold (from policy by default).
            max_boundary_slack: Maximum terminal slack threshold (from policy by default).
            min_mapq: Minimum mapping quality threshold (from policy by default).
            reads_sequences: Optional dictionary mapping read ID to nucleotide sequence.

        Returns:
            List of validated SpanningRead objects sorted by read ID.
        """
        anchor_min = min_anchor_bp if min_anchor_bp is not None else self.policy.ont.min_anchor_bp.value
        slack_max = (
            max_boundary_slack if max_boundary_slack is not None else self.policy.ont.max_boundary_slack_bp.value
        )
        mapq_min = min_mapq if min_mapq is not None else self.policy.ont.min_mapq.value

        left_name, left_strand = _parse_port_spec(left_port)
        right_name, right_strand = _parse_port_spec(right_port)

        left_sym = left_port.to_symbol() if hasattr(left_port, "to_symbol") else str(left_port)
        right_sym = right_port.to_symbol() if hasattr(right_port, "to_symbol") else str(right_port)

        left_match_names = {left_name, left_sym, f"{left_name}+", f"{left_name}-"}
        right_match_names = {right_name, right_sym, f"{right_name}+", f"{right_name}-"}
        target_names = left_match_names | right_match_names

        # Memory-efficient collection: only retain hits matching target ports
        hits_by_read: dict[str, list[PafHit]] = defaultdict(list)
        source_iterable: Iterable[str | PafHit]
        if isinstance(paf_source, str | Path | TextIO):
            source_iterable = parse_paf_hits(paf_source)
        else:
            source_iterable = paf_source

        for item in source_iterable:
            if isinstance(item, PafHit):
                hit = item
            else:
                s = item.strip()
                if not s or s.startswith("#"):
                    continue
                hit = PafHit.from_line(s)
            if hit.tname in target_names:
                hits_by_read[hit.qname].append(hit)

        spanning_reads: list[SpanningRead] = []

        fwd_l_strand = left_strand.to_symbol()
        fwd_r_strand = right_strand.to_symbol()
        rev_l_strand = left_strand.invert().to_symbol()
        rev_r_strand = right_strand.invert().to_symbol()

        for qname, read_hits in hits_by_read.items():
            l_hits = [
                h
                for h in read_hits
                if h.tname in left_match_names
                and h.alen >= anchor_min
                and h.mapq >= mapq_min
                and _port_boundary_slack(h, is_left=True, strand=left_strand) <= slack_max
            ]
            r_hits = [
                h
                for h in read_hits
                if h.tname in right_match_names
                and h.alen >= anchor_min
                and h.mapq >= mapq_min
                and _port_boundary_slack(h, is_left=False, strand=right_strand) <= slack_max
            ]

            if not l_hits or not r_hits:
                continue

            best_candidate: SpanningRead | None = None
            min_inner_distance: int | None = None

            # 1. Forward walk orientation: Left before Right (L.qend <= R.qstart)
            l_fwd = [h for h in l_hits if h.strand == fwd_l_strand]
            r_fwd = [h for h in r_hits if h.strand == fwd_r_strand]
            for lh in l_fwd:
                for rh in r_fwd:
                    if lh.qend <= rh.qstart:
                        l_slack = _port_boundary_slack(lh, is_left=True, strand=left_strand)
                        r_slack = _port_boundary_slack(rh, is_left=False, strand=right_strand)
                        begin = lh.qend + l_slack
                        end = rh.qstart - r_slack
                        if end >= begin:
                            inner_dist = rh.qstart - lh.qend
                            if min_inner_distance is None or inner_dist < min_inner_distance:
                                min_inner_distance = inner_dist
                                seq = None
                                if reads_sequences and qname in reads_sequences:
                                    seq = reads_sequences[qname][begin:end]
                                best_candidate = SpanningRead(
                                    qname=qname,
                                    strand="+",
                                    left_hit=lh,
                                    right_hit=rh,
                                    gap_length=end - begin,
                                    q_from=begin,
                                    q_to=end,
                                    sequence=seq,
                                )

            # 2. Reverse walk orientation: Right before Left (R.qend <= L.qstart)
            l_rev = [h for h in l_hits if h.strand == rev_l_strand]
            r_rev = [h for h in r_hits if h.strand == rev_r_strand]
            for lh in l_rev:
                for rh in r_rev:
                    if rh.qend <= lh.qstart:
                        r_slack = _port_boundary_slack(rh, is_left=False, strand=right_strand)
                        l_slack = _port_boundary_slack(lh, is_left=True, strand=left_strand)
                        begin = rh.qend + r_slack
                        end = lh.qstart - l_slack
                        if end >= begin:
                            inner_dist = lh.qstart - rh.qend
                            if min_inner_distance is None or inner_dist < min_inner_distance:
                                min_inner_distance = inner_dist
                                seq = None
                                if reads_sequences and qname in reads_sequences:
                                    raw_seq = reads_sequences[qname][begin:end]
                                    seq = reverse_complement(raw_seq)
                                best_candidate = SpanningRead(
                                    qname=qname,
                                    strand="-",
                                    left_hit=lh,
                                    right_hit=rh,
                                    gap_length=end - begin,
                                    q_from=begin,
                                    q_to=end,
                                    sequence=seq,
                                )

            if best_candidate is not None:
                spanning_reads.append(best_candidate)

        spanning_reads.sort(key=lambda r: r.qname)
        return spanning_reads

    def compute_distribution(
        self,
        spanning_reads_or_lengths: Sequence[SpanningRead] | Sequence[int],
    ) -> GapDistribution:
        """Compute empirical minimum, median, maximum, and IQR across spanning reads.

        Args:
            spanning_reads_or_lengths: Sequence of SpanningRead objects or raw integer lengths.

        Returns:
            GapDistribution with computed summary statistics.
        """
        if not spanning_reads_or_lengths:
            return GapDistribution(min_bp=0, median_bp=0, max_bp=0, iqr_bp=0, count=0)

        lengths: list[int] = []
        for item in spanning_reads_or_lengths:
            if isinstance(item, SpanningRead):
                lengths.append(item.gap_length)
            elif isinstance(item, int):
                lengths.append(item)
            else:
                lengths.append(int(getattr(item, "gap_length", item)))

        if not lengths:
            return GapDistribution(min_bp=0, median_bp=0, max_bp=0, iqr_bp=0, count=0)

        lengths.sort()
        count = len(lengths)
        min_bp = lengths[0]
        max_bp = lengths[-1]
        med = statistics.median(lengths)
        median_bp = round(med)

        if count >= 2:
            n_q = self.policy.ont.quantile_num.value
            qs = statistics.quantiles(lengths, n=n_q)
            q1 = qs[0]
            q3 = qs[-1]
            iqr_bp = round(q3 - q1)
        else:
            iqr_bp = 0

        return GapDistribution(
            min_bp=min_bp,
            median_bp=median_bp,
            max_bp=max_bp,
            iqr_bp=iqr_bp,
            count=count,
        )

    def resolve_satellite_gap(
        self,
        chrom: str,
        side: str,
        left_port: str | SegmentPort,
        right_port: str | SegmentPort,
        paf_source: Path | str | TextIO | Iterable[str] | Iterable[PafHit] | None = None,
        *,
        gap_index: int = 1,
        spanning_reads: Sequence[SpanningRead] | None = None,
        reads_sequences: Mapping[str, str] | None = None,
        repeat_monomer_len_bp: int | None = None,
        core_coverage_bp_median: float | None = None,
        core_depth_ratio: float | None = None,
        somatic_depth_ratio: float | None = None,
        somatic_port: str | None = None,
        matched_end: EndStateRecord | None = None,
        min_spanning_reads: int | None = None,
        notes: str = "",
    ) -> SatelliteGapResolution:
        """Resolve a variable satellite gap locus to END_COMPLETE or END_GERMLINE_ORPHAN.

        Implements BIO-004:
        - When bridged by >= min_spanning_reads, selects representative sequence from
          candidate read closest to empirical median gap length without generic N collapse.
        - Resolves chromosome end state to EndState.END_COMPLETE (or END_GERMLINE_ORPHAN if unbridged).
        - Emits typed QuirkRecord and EndStateRecord instances for biological ledger tracking.

        Args:
            chrom: Chromosome name (e.g. 'IV', 'V').
            side: Arm end ('5p' or '3p').
            left_port: Proximal unitig port bounding gap.
            right_port: Distal unitig port bounding gap.
            paf_source: Optional PAF file path or stream to extract spanning reads.
            gap_index: 1-based occurrence index for the quirk feature identifier.
            spanning_reads: Pre-collected SpanningRead sequence if already computed.
            reads_sequences: Optional dictionary mapping read ID to sequence.
            repeat_monomer_len_bp: Known repeating satellite monomer length in bp.
            core_coverage_bp_median: Somatic core reference depth ($C_{core}$).
            core_depth_ratio: Relative depth ratio ($D / C_{core}$).
            somatic_depth_ratio: Optional somatic cap depth ratio relative to core.
            min_spanning_reads: Minimum spanning reads threshold (from policy by default).
            notes: Contextual provenance notes.

        Returns:
            SatelliteGapResolution containing typed QuirkRecord, EndStateRecord, and candidates.
        """
        reads_list: list[SpanningRead]
        if spanning_reads is None:
            if paf_source is None:
                reads_list = []
            else:
                reads_list = self.collect_spanning_reads(
                    paf_source,
                    left_port,
                    right_port,
                    reads_sequences=reads_sequences,
                )
        else:
            reads_list = list(spanning_reads)

        floor = min_spanning_reads if min_spanning_reads is not None else self.policy.ont.min_spanning_reads.value

        left_str = left_port.to_symbol() if hasattr(left_port, "to_symbol") else str(left_port)
        right_str = right_port.to_symbol() if hasattr(right_port, "to_symbol") else str(right_port)
        left_name = getattr(left_port, "name", str(left_port))
        right_name = getattr(right_port, "name", str(right_port))

        actual_somatic_port = left_str
        if matched_end is not None and matched_end.somatic_port is not None:
            actual_somatic_port = matched_end.somatic_port
        elif somatic_port is not None:
            actual_somatic_port = somatic_port

        feature_id_quirk = build_feature_id(IntervalQuirk.SAT_GAP_VARIABLE, chrom, side, index=gap_index)

        if len(reads_list) >= floor:
            # Successfully bridged
            dist = self.compute_distribution(reads_list)
            best_read = min(
                reads_list,
                key=lambda c: (abs(c.gap_length - dist.median_bp), c.gap_length, c.qname),
            )

            rep_seq = best_read.sequence
            if rep_seq is None and reads_sequences and best_read.qname in reads_sequences:
                raw_seq = reads_sequences[best_read.qname][best_read.q_from : best_read.q_to]
                rep_seq = reverse_complement(raw_seq) if best_read.strand == "-" else raw_seq

            if rep_seq is not None and best_read.sequence is None:
                best_read = SpanningRead(
                    qname=best_read.qname,
                    strand=best_read.strand,
                    left_hit=best_read.left_hit,
                    right_hit=best_read.right_hit,
                    gap_length=best_read.gap_length,
                    q_from=best_read.q_from,
                    q_to=best_read.q_to,
                    sequence=rep_seq,
                )

            copy_number_range: str | None = None
            if repeat_monomer_len_bp and repeat_monomer_len_bp > 0:
                cn_min = dist.min_bp // repeat_monomer_len_bp
                cn_max = dist.max_bp // repeat_monomer_len_bp
                copy_number_range = f"{cn_min}-{cn_max}"

            state = EndState.END_COMPLETE
            feature_id_end = build_feature_id(state, chrom, side)

            quirk_notes = (
                f"Bridged by {len(reads_list)} ONT reads; representative sequence selected "
                f"from {best_read.qname} closest to median gap {dist.median_bp} bp (BIO-004)"
            )
            if notes:
                quirk_notes = f"{quirk_notes}; {notes}"

            end_notes = (
                f"Resolved via ONT bridging across variable satellite gap with {len(reads_list)} "
                f"spanning reads; median gap {dist.median_bp} bp (BIO-004)"
            )
            if notes:
                end_notes = f"{end_notes}; {notes}"

            quirk_record = QuirkRecord(
                feature_id=feature_id_quirk,
                chrom=chrom,
                side=side,
                quirk_type=IntervalQuirk.SAT_GAP_VARIABLE,
                spanning_read_count=len(reads_list),
                length_min_bp=dist.min_bp,
                length_median_bp=dist.median_bp,
                length_max_bp=dist.max_bp,
                length_iqr_bp=dist.iqr_bp,
                repeat_monomer_len_bp=repeat_monomer_len_bp,
                copy_number_range=copy_number_range,
                left_anchor_port=left_str,
                right_anchor_port=right_str,
                agp_component_type="W",
                agp_linkage_evidence="ont-spanning",
                unitig_ids=(left_name, right_name),
                read_count=len(reads_list),
                core_coverage_bp_median=core_coverage_bp_median,
                core_depth_ratio=core_depth_ratio,
                representative_sequence_selected=True,
                sequence_provenance="ONT-derived",
                unpolished=True,
                notes=quirk_notes,
            )

            end_state_record = EndStateRecord(
                feature_id=feature_id_end,
                chrom=chrom,
                side=side,
                state=state,
                somatic_port=actual_somatic_port,
                germline_port=right_str,
                shared_port=None,
                core_coverage_bp_median=core_coverage_bp_median,
                somatic_depth_ratio=somatic_depth_ratio,
                germline_depth_ratio=core_depth_ratio,
                resolved=True,
                notes=end_notes,
            )

            res = SatelliteGapResolution(
                chrom=chrom,
                side=side,
                state=state,
                quirk_record=quirk_record,
                end_state_record=end_state_record,
                spanning_reads=reads_list,
                representative_read=best_read,
                representative_sequence=rep_seq,
                distribution=dist,
            )
            self.gap_resolutions.append(res)
            return res

        # Unbridged orphan cap
        state = EndState.END_GERMLINE_ORPHAN
        feature_id_end = build_feature_id(state, chrom, side)
        dist = self.compute_distribution(reads_list) if reads_list else None

        quirk_notes = f"Unbridged variable satellite gap: {len(reads_list)} spanning reads found (policy floor {floor})"
        if notes:
            quirk_notes = f"{quirk_notes}; {notes}"

        end_notes = f"Unbridged germline orphan cap at {chrom}:{side}; insufficient ONT spanning support"
        if notes:
            end_notes = f"{end_notes}; {notes}"

        quirk_record = QuirkRecord(
            feature_id=feature_id_quirk,
            chrom=chrom,
            side=side,
            quirk_type=IntervalQuirk.SAT_GAP_VARIABLE,
            spanning_read_count=len(reads_list),
            length_min_bp=dist.min_bp if dist else None,
            length_median_bp=dist.median_bp if dist else None,
            length_max_bp=dist.max_bp if dist else None,
            length_iqr_bp=dist.iqr_bp if dist else None,
            repeat_monomer_len_bp=repeat_monomer_len_bp,
            copy_number_range=None,
            left_anchor_port=left_str,
            right_anchor_port=right_str,
            agp_component_type="U",
            agp_linkage_evidence="unbridged",
            unitig_ids=(left_name, right_name),
            read_count=len(reads_list),
            core_coverage_bp_median=core_coverage_bp_median,
            core_depth_ratio=core_depth_ratio,
            representative_sequence_selected=False,
            notes=quirk_notes,
        )

        end_state_record = EndStateRecord(
            feature_id=feature_id_end,
            chrom=chrom,
            side=side,
            state=state,
            somatic_port=actual_somatic_port,
            germline_port=None,
            shared_port=None,
            core_coverage_bp_median=core_coverage_bp_median,
            somatic_depth_ratio=somatic_depth_ratio,
            germline_depth_ratio=core_depth_ratio,
            resolved=False,
            notes=end_notes,
        )

        res = SatelliteGapResolution(
            chrom=chrom,
            side=side,
            state=state,
            quirk_record=quirk_record,
            end_state_record=end_state_record,
            spanning_reads=reads_list,
            representative_read=None,
            representative_sequence=None,
            distribution=dist,
        )
        self.gap_resolutions.append(res)
        return res

    def _find_distal_germline_targets(self, germline_port: SegmentPort) -> list[SegmentPort]:
        """Find candidate distal germline telomeric caps reachable from germline_port."""
        if self.graph is None:
            return []
        distal_targets: list[SegmentPort] = []
        visited = {germline_port.name}
        queue = [germline_port]
        while queue:
            p = queue.pop(0)
            exit_port = p.invert()
            outgoing = self.graph.get_outgoing(exit_port)
            for edge in outgoing:
                nxt = edge.to_port
                if nxt.name not in visited:
                    visited.add(nxt.name)
                    distal_targets.append(nxt)
                    if self.graph.degree(nxt.invert()) > 0:
                        queue.append(nxt)
        orphan_caps = self.identify_orphan_caps()
        for cap_seg, (_tel_port, inward_port) in orphan_caps.items():
            if cap_seg not in visited:
                p_name, p_strand = _parse_port_spec(inward_port)
                distal_targets.append(SegmentPort(name=p_name, strand=p_strand))
        return distal_targets

    def resolve_somatic_dropout(
        self,
        chrom: str,
        side: str,
        stem_port: SegmentPort,
        germline_port: SegmentPort,
        somatic_port: SegmentPort,
        paf_path: str | Path,
        reads_fasta: dict[str, str] | None = None,
        *,
        core_coverage_bp_median: float | None = None,
        somatic_depth_ratio: float | None = None,
        germline_depth_ratio: float | None = None,
        core_depth_ratio: float | None = None,
        gap_index: int = 1,
        notes: str = "",
    ) -> tuple[QuirkRecord, EndStateRecord, str | None]:
        """Resolve a somatic dropout chromosome end via spanning ONT reads.

        At chromosome ends where the somatic branch suffers a sequencing dropout,
        searches for spanning ONT reads bridging the stem port into the germline
        cap unitig (or through to the distal germline telomeric cap). If spanning
        reads exist:
        - Extracts representative sequence from candidate read closest to empirical
          median gap length (BIO-004).
        - Resolves chromosome end state to EndState.END_COMPLETE (resolved=True),
          setting notes="ONT-derived sequence resolved germline telomere".
        - Emits QuirkRecord documenting the rescued dropout stub
          (IntervalQuirk.LOW_SUPPORT_STUB), with representative_sequence_selected=True.

        If no spanning reads exist:
        - Retains EndState.END_SOMATIC_MISSING (resolved=False).
        - Emits QuirkRecord with representative_sequence_selected=False.

        Args:
            chrom: Chromosome name (e.g. 'IV', 'V').
            side: Arm end ('5p' or '3p').
            stem_port: Stem unitig port entering the bifurcation.
            germline_port: Intact germline cap or continuation unitig port.
            somatic_port: Somatic dropout stub unitig port.
            paf_path: PAF file path or stream to extract spanning reads.
            reads_fasta: Optional dictionary mapping read ID to nucleotide sequence.
            core_coverage_bp_median: Somatic core reference depth ($C_{core}$).
            somatic_depth_ratio: Somatic cap depth ratio relative to core.
            germline_depth_ratio: Germline cap depth ratio relative to core.
            core_depth_ratio: Feature depth relative to somatic core coverage.
            gap_index: 1-based occurrence index for quirk feature identifier.
            notes: Contextual provenance notes.

        Returns:
            Tuple of (QuirkRecord, EndStateRecord, representative_sequence).
        """
        cand_targets = [germline_port]
        cand_targets.extend(self._find_distal_germline_targets(germline_port))

        spanning_reads: list[SpanningRead] = []
        resolved_germline_port = germline_port

        if isinstance(paf_path, str | Path):
            p = Path(paf_path)
            if p.exists():
                for target in cand_targets:
                    reads = self.collect_spanning_reads(
                        p,
                        stem_port,
                        target,
                        reads_sequences=reads_fasta,
                    )
                    if reads:
                        spanning_reads = reads
                        resolved_germline_port = target
                        break
        else:
            for target in cand_targets:
                reads = self.collect_spanning_reads(
                    paf_path,
                    stem_port,
                    target,
                    reads_sequences=reads_fasta,
                )
                if reads:
                    spanning_reads = reads
                    resolved_germline_port = target
                    break

        stem_str = stem_port.to_symbol() if hasattr(stem_port, "to_symbol") else str(stem_port)
        germline_str = (
            resolved_germline_port.to_symbol()
            if hasattr(resolved_germline_port, "to_symbol")
            else str(resolved_germline_port)
        )
        somatic_str = somatic_port.to_symbol() if hasattr(somatic_port, "to_symbol") else str(somatic_port)
        somatic_name = getattr(somatic_port, "name", str(somatic_port))

        feature_id_quirk = build_feature_id(IntervalQuirk.LOW_SUPPORT_STUB, chrom, side, index=gap_index)
        core_ratio = core_depth_ratio if core_depth_ratio is not None else somatic_depth_ratio

        if spanning_reads:
            dist = self.compute_distribution(spanning_reads)
            best_read = min(
                spanning_reads,
                key=lambda c: (abs(c.gap_length - dist.median_bp), c.gap_length, c.qname),
            )

            rep_seq = best_read.sequence
            if rep_seq is None and reads_fasta and best_read.qname in reads_fasta:
                raw_seq = reads_fasta[best_read.qname][best_read.q_from : best_read.q_to]
                rep_seq = reverse_complement(raw_seq) if best_read.strand == "-" else raw_seq

            state = EndState.END_COMPLETE
            feature_id_end = build_feature_id(state, chrom, side)

            end_notes = "ONT-derived sequence resolved germline telomere"
            if notes:
                end_notes = f"{end_notes}; {notes}"

            quirk_notes = (
                f"Rescued somatic dropout stub via {len(spanning_reads)} spanning ONT reads; "
                f"unpolished ONT-derived representative sequence selected from {best_read.qname} "
                f"closest to median gap {dist.median_bp} bp (BIO-004)"
            )
            if notes:
                quirk_notes = f"{quirk_notes}; {notes}"

            quirk_record = QuirkRecord(
                feature_id=feature_id_quirk,
                chrom=chrom,
                side=side,
                quirk_type=IntervalQuirk.LOW_SUPPORT_STUB,
                spanning_read_count=len(spanning_reads),
                length_min_bp=dist.min_bp,
                length_median_bp=dist.median_bp,
                length_max_bp=dist.max_bp,
                length_iqr_bp=dist.iqr_bp,
                left_anchor_port=stem_str,
                right_anchor_port=germline_str,
                agp_component_type="W",
                agp_linkage_evidence="ont-spanning",
                unitig_ids=(somatic_name,),
                read_count=len(spanning_reads),
                core_coverage_bp_median=core_coverage_bp_median,
                core_depth_ratio=core_ratio,
                representative_sequence_selected=True,
                sequence_provenance="ONT-derived",
                unpolished=True,
                notes=quirk_notes,
            )

            end_state_record = EndStateRecord(
                feature_id=feature_id_end,
                chrom=chrom,
                side=side,
                state=state,
                somatic_port=somatic_str,
                germline_port=germline_str,
                shared_port=None,
                core_coverage_bp_median=core_coverage_bp_median,
                somatic_depth_ratio=somatic_depth_ratio,
                germline_depth_ratio=germline_depth_ratio,
                resolved=True,
                notes=end_notes,
            )

            res_dropout = (quirk_record, end_state_record, rep_seq)
            self.dropout_resolutions.append(res_dropout)
            return res_dropout

        state = EndState.END_SOMATIC_MISSING
        feature_id_end = build_feature_id(state, chrom, side)

        end_notes = f"Unrescued somatic dropout at {chrom}:{side}; insufficient ONT spanning support"
        if notes:
            end_notes = f"{end_notes}; {notes}"

        quirk_notes = f"Unrescued somatic dropout stub at {chrom}:{side}; no spanning ONT reads found"
        if notes:
            quirk_notes = f"{quirk_notes}; {notes}"

        quirk_record = QuirkRecord(
            feature_id=feature_id_quirk,
            chrom=chrom,
            side=side,
            quirk_type=IntervalQuirk.LOW_SUPPORT_STUB,
            spanning_read_count=0,
            length_min_bp=None,
            length_median_bp=None,
            length_max_bp=None,
            length_iqr_bp=None,
            left_anchor_port=stem_str,
            right_anchor_port=germline_str,
            agp_component_type="U",
            agp_linkage_evidence="unbridged",
            unitig_ids=(somatic_name,),
            read_count=0,
            core_coverage_bp_median=core_coverage_bp_median,
            core_depth_ratio=core_ratio,
            representative_sequence_selected=False,
            notes=quirk_notes,
        )

        end_state_record = EndStateRecord(
            feature_id=feature_id_end,
            chrom=chrom,
            side=side,
            state=state,
            somatic_port=None,
            germline_port=germline_str,
            shared_port=None,
            core_coverage_bp_median=core_coverage_bp_median,
            somatic_depth_ratio=somatic_depth_ratio,
            germline_depth_ratio=germline_depth_ratio,
            resolved=False,
            notes=end_notes,
        )

        res_unrescued = (quirk_record, end_state_record, None)
        self.dropout_resolutions.append(res_unrescued)
        return res_unrescued
