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
from collections import defaultdict
from collections.abc import Iterable, Iterator, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import TextIO

from gfapdefinisher.models import SegmentPort
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


class OntBridgingEngine:
    """Resolution engine bridging assembly graph ports with long ONT reads.

    Provides innermost boundary hit detection to resolve multivalent tandem satellite
    repeat alignments, empirical gap length distribution statistics, and BIO-004
    representative sequence selection.
    """

    def __init__(self, policy: Policy = DEFAULT_POLICY) -> None:
        self.policy = policy

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

        left_name = left_port.name if isinstance(left_port, SegmentPort) else str(left_port)
        right_name = right_port.name if isinstance(right_port, SegmentPort) else str(right_port)
        left_sym = left_port.to_symbol() if hasattr(left_port, "to_symbol") else str(left_port)
        right_sym = right_port.to_symbol() if hasattr(right_port, "to_symbol") else str(right_port)

        target_names = {left_name, right_name, left_sym, right_sym}

        # Memory-efficient collection: only retain hits matching target ports
        hits_by_read: dict[str, list[PafHit]] = defaultdict(list)
        if isinstance(paf_source, Iterable) and not isinstance(paf_source, str | Path | TextIO):
            # Could be PafHit stream or raw line stream
            first_item = next(iter(paf_source), None)
            if first_item is None:
                return []
            if isinstance(first_item, PafHit):
                if first_item.tname in target_names:
                    hits_by_read[first_item.qname].append(first_item)
                for h in paf_source:
                    if isinstance(h, PafHit) and h.tname in target_names:
                        hits_by_read[h.qname].append(h)
            else:
                for line in paf_source:
                    s = line.strip()
                    if s and not s.startswith("#"):
                        hit = PafHit.from_line(s)
                        if hit.tname in target_names:
                            hits_by_read[hit.qname].append(hit)
        else:
            for hit in parse_paf_hits(paf_source):
                if hit.tname in target_names:
                    hits_by_read[hit.qname].append(hit)

        spanning_reads: list[SpanningRead] = []

        for qname, read_hits in hits_by_read.items():
            l_hits = [
                h
                for h in read_hits
                if h.tname in (left_name, left_sym)
                and h.alen >= anchor_min
                and h.mapq >= mapq_min
                and h.left_slack <= slack_max
            ]
            r_hits = [
                h
                for h in read_hits
                if h.tname in (right_name, right_sym)
                and h.alen >= anchor_min
                and h.mapq >= mapq_min
                and h.right_slack <= slack_max
            ]

            if not l_hits or not r_hits:
                continue

            best_candidate: SpanningRead | None = None
            min_inner_distance: int | None = None

            # 1. Forward strand dovetail collinearity: Left before Right (L.qend <= R.qstart)
            l_fwd = [h for h in l_hits if h.strand == "+"]
            r_fwd = [h for h in r_hits if h.strand == "+"]
            for lh in l_fwd:
                for rh in r_fwd:
                    if lh.qend <= rh.qstart:
                        begin = lh.qend + lh.left_slack
                        end = rh.qstart - rh.right_slack
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

            # 2. Reverse strand dovetail collinearity: Right before Left (R.qend <= L.qstart)
            l_rev = [h for h in l_hits if h.strand == "-"]
            r_rev = [h for h in r_hits if h.strand == "-"]
            for lh in l_rev:
                for rh in r_rev:
                    if rh.qend <= lh.qstart:
                        begin = rh.qend + rh.right_slack
                        end = lh.qstart - lh.left_slack
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
        spanning_reads: Sequence[SpanningRead] | None = None,
        reads_sequences: Mapping[str, str] | None = None,
        repeat_monomer_len_bp: int | None = None,
        core_coverage_bp_median: float | None = None,
        core_depth_ratio: float | None = None,
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
            spanning_reads: Pre-collected SpanningRead sequence if already computed.
            reads_sequences: Optional dictionary mapping read ID to sequence.
            repeat_monomer_len_bp: Known repeating satellite monomer length in bp.
            core_coverage_bp_median: Somatic core reference depth ($C_{core}$).
            core_depth_ratio: Relative depth ratio ($D / C_{core}$).
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

        feature_id_quirk = build_feature_id(IntervalQuirk.SAT_GAP_VARIABLE, chrom, side, index=1)

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
                notes=quirk_notes,
            )

            end_state_record = EndStateRecord(
                feature_id=feature_id_end,
                chrom=chrom,
                side=side,
                state=state,
                somatic_port=left_str,
                germline_port=right_str,
                shared_port=None,
                core_coverage_bp_median=core_coverage_bp_median,
                germline_depth_ratio=core_depth_ratio,
                resolved=True,
                notes=end_notes,
            )

            return SatelliteGapResolution(
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
            somatic_port=left_str,
            germline_port=None,
            shared_port=None,
            core_coverage_bp_median=core_coverage_bp_median,
            germline_depth_ratio=core_depth_ratio,
            resolved=False,
            notes=end_notes,
        )

        return SatelliteGapResolution(
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
