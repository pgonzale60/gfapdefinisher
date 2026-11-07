"""Precomputed annotation layer loaders for assembly graph nodes.

Provides data structures and loaders for external evidence layers:
- Tiara: deep-learning taxonomic / contaminant classification
- TIDK: telomeric repeat window counts
- Barrnap: ribosomal RNA feature predictions (5S, 5.8S, 18S, 28S)
- Compleasm: lineage-specific USCO single-copy ortholog hits
- NodeAnnotations: integrated multi-layer container projected per unitig
"""

from __future__ import annotations

import csv
import gzip
from collections.abc import Iterable
from dataclasses import dataclass, field
from pathlib import Path


@dataclass(slots=True, frozen=True)
class TiaraHit:
    """Taxonomic or contaminant classification hit from Tiara."""

    sequence_id: str
    classification: str


@dataclass(slots=True, frozen=True)
class TidkWindow:
    """Telomeric repeat motif count within a sliding sequence window from TIDK.

    Note: TIDK reports the 1-based end coordinate of each window in the 'window' column.
    """

    sequence_id: str
    window_end: int
    forward_count: int
    reverse_count: int
    window_start: int

    def __init__(
        self,
        sequence_id: str,
        window_end: int | None = None,
        forward_count: int = 0,
        reverse_count: int = 0,
        window_start: int | None = None,
    ) -> None:
        if window_end is None and window_start is not None:
            end_val = window_start
            start_val = window_start
        elif window_end is not None and window_start is None:
            end_val = window_end
            start_val = window_end
        elif window_end is not None and window_start is not None:
            end_val = window_end
            start_val = window_start
        else:
            end_val = 0
            start_val = 0
        object.__setattr__(self, "sequence_id", sequence_id)
        object.__setattr__(self, "window_end", end_val)
        object.__setattr__(self, "forward_count", forward_count)
        object.__setattr__(self, "reverse_count", reverse_count)
        object.__setattr__(self, "window_start", start_val)


TidkWindowRow = TidkWindow


@dataclass(slots=True, frozen=True)
class BarrnapHit:
    """Ribosomal RNA gene feature prediction from Barrnap."""

    sequence_id: str
    start: int
    end: int
    strand: str
    feature_name: str
    evalue: float

    @property
    def product(self) -> str:
        """Alias for feature_name matching GFF product attribute."""
        return self.feature_name


@dataclass(slots=True, frozen=True)
class CompleasmHit:
    """Universal single-copy ortholog (USCO) benchmark hit from Compleasm."""

    gene_id: str
    status: str
    sequence_id: str
    start: int
    end: int
    score: float


@dataclass(slots=True)
class NodeAnnotations:
    """Container holding multi-layer precomputed annotations for an assembly node."""

    sequence_id: str
    tiara: TiaraHit | None = None
    tidk_windows: list[TidkWindow] = field(default_factory=list)
    barrnap_hits: list[BarrnapHit] = field(default_factory=list)
    compleasm_hits: list[CompleasmHit] = field(default_factory=list)

    def __post_init__(self) -> None:
        if isinstance(self.tiara, str):
            self.tiara = TiaraHit(sequence_id=self.sequence_id, classification=self.tiara)

    @property
    def tiara_classification(self) -> str | None:
        """Return Tiara classification string if available."""
        return self.tiara.classification if self.tiara is not None else None

    @property
    def tiara_hit(self) -> TiaraHit | None:
        """Alias for tiara attribute."""
        return self.tiara

    @property
    def tidk(self) -> list[TidkWindow]:
        """Alias for tidk_windows."""
        return self.tidk_windows

    @property
    def barrnap(self) -> list[BarrnapHit]:
        """Alias for barrnap_hits."""
        return self.barrnap_hits

    @property
    def compleasm(self) -> list[CompleasmHit]:
        """Alias for compleasm_hits."""
        return self.compleasm_hits

    def has_telomere(self, threshold: int = 1) -> bool:
        """Return True if any TIDK window contains repeat counts >= threshold."""
        return any(w.forward_count >= threshold or w.reverse_count >= threshold for w in self.tidk_windows)

    def has_rdna(self) -> bool:
        """Return True if node has one or more Barrnap rDNA annotations."""
        return bool(self.barrnap_hits)

    def is_eukaryotic(self) -> bool:
        """Return True if Tiara classification identifies node as eukaryotic."""
        return self.tiara is not None and self.tiara.classification.lower() == "eukarya"

    def is_contaminant(self) -> bool:
        """Return True if Tiara classification identifies node as prokaryotic/contaminant."""
        return self.tiara is not None and self.tiara.classification.lower() in (
            "bacteria",
            "archaea",
            "prokarya",
            "contaminant",
        )


def _open_file(path: str | Path):
    """Open a text file path supporting optional gzip compression."""
    p = Path(path)
    if p.suffix == ".gz":
        return gzip.open(p, mode="rt", encoding="utf-8")
    return open(p, encoding="utf-8")


def load_tiara(path: str | Path) -> dict[str, str]:
    """Parse Tiara classification table (tiara_unitigs.txt).

    Args:
        path: Path to Tiara tabular output.

    Returns:
        Mapping from sequence_id to classification label (e.g. 'eukarya', 'bacteria').
    """
    classifications: dict[str, str] = {}
    with _open_file(path) as f:
        reader = csv.reader(f, delimiter="\t")
        header = None
        seq_idx = -1
        class_idx = -1

        for row in reader:
            if not row:
                continue
            first_val = row[0].strip()
            if not first_val or first_val.startswith("#"):
                continue

            if header is None:
                header = [col.strip().lower() for col in row]
                for idx, col in enumerate(header):
                    if col in ("sequence_id", "id", "seq_id", "sequence"):
                        seq_idx = idx
                    elif col in ("class_fst_stage", "classification", "class"):
                        class_idx = idx
                if seq_idx == -1:
                    seq_idx = 0
                if class_idx == -1:
                    class_idx = 1
                continue

            if len(row) <= max(seq_idx, class_idx):
                continue

            seq_id = row[seq_idx].strip()
            classification = row[class_idx].strip()
            if seq_id:
                classifications[seq_id] = classification

    return classifications


def load_tidk(path: str | Path) -> dict[str, list[TidkWindow]]:
    """Parse TIDK telomeric repeat windows TSV.

    Args:
        path: Path to tidk_telomeres_telomeric_repeat_windows.tsv.

    Returns:
        Mapping from sequence_id to list of TidkWindow records in window order.
    """
    windows: dict[str, list[TidkWindow]] = {}
    with _open_file(path) as f:
        reader = csv.reader(f, delimiter="\t")
        header = None
        seq_idx = -1
        win_idx = -1
        fwd_idx = -1
        rev_idx = -1

        for row in reader:
            if not row:
                continue
            first_val = row[0].strip()
            if not first_val or first_val.startswith("#"):
                continue

            if header is None:
                header = [col.strip().lower() for col in row]
                for idx, col in enumerate(header):
                    if col in ("id", "sequence_id", "seq_id", "sequence"):
                        seq_idx = idx
                    elif col in ("window", "window_end", "window_start", "start", "end"):
                        win_idx = idx
                    elif col in ("forward_repeat_number", "forward_count", "forward", "fwd"):
                        fwd_idx = idx
                    elif col in ("reverse_repeat_number", "reverse_count", "reverse", "rev"):
                        rev_idx = idx

                if seq_idx == -1:
                    seq_idx = 0
                if win_idx == -1:
                    win_idx = 1
                if fwd_idx == -1:
                    fwd_idx = 2
                if rev_idx == -1:
                    rev_idx = 3
                continue

            max_idx = max(seq_idx, win_idx, fwd_idx, rev_idx)
            if len(row) <= max_idx:
                continue

            seq_id = row[seq_idx].strip()
            if not seq_id:
                continue

            hit = TidkWindow(
                sequence_id=seq_id,
                window_end=int(float(row[win_idx].strip())),
                forward_count=int(float(row[fwd_idx].strip())),
                reverse_count=int(float(row[rev_idx].strip())),
            )
            windows.setdefault(seq_id, []).append(hit)

    return windows


def load_barrnap(path: str | Path) -> dict[str, list[BarrnapHit]]:
    """Parse Barrnap rDNA annotations in GFF3 format.

    Args:
        path: Path to barrnap_rdna.gff.

    Returns:
        Mapping from sequence_id to list of BarrnapHit records.
    """
    gff_fields = ("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
    hits: dict[str, list[BarrnapHit]] = {}

    with _open_file(path) as f:
        for line in f:
            line_str = line.strip()
            if not line_str or line_str.startswith("#"):
                continue

            parts = [p.strip() for p in line_str.split("\t")]
            while len(parts) < len(gff_fields):
                parts.append("")

            seqid, _src, type_name, start_s, end_s, score_s, strand_s, _phase, attr_s = parts[: len(gff_fields)]
            if not seqid:
                continue

            start = int(float(start_s)) if start_s else 0
            end = int(float(end_s)) if end_s else 0
            evalue = float(score_s) if score_s not in (".", "") else float(0)
            strand = strand_s if strand_s in ("+", "-") else "+"

            feature_name = ""
            if attr_s:
                attrs: dict[str, str] = {}
                for item in attr_s.split(";"):
                    if "=" in item:
                        k, v = item.split("=", 1)
                        attrs[k.strip()] = v.strip()
                feature_name = attrs.get("Name") or attrs.get("ID") or attrs.get("product") or ""

            if not feature_name:
                feature_name = type_name or "rRNA"

            hit = BarrnapHit(
                sequence_id=seqid,
                start=start,
                end=end,
                strand=strand,
                feature_name=feature_name,
                evalue=evalue,
            )
            hits.setdefault(seqid, []).append(hit)

    return hits


def load_compleasm(path: str | Path) -> dict[str, list[CompleasmHit]]:
    """Parse Compleasm USCO full table TSV (compleasm_unitigs_full_table.tsv).

    Args:
        path: Path to Compleasm full table TSV.

    Returns:
        Mapping from sequence_id to list of CompleasmHit records placed on that unitig.
    """

    hits: dict[str, list[CompleasmHit]] = {}
    with _open_file(path) as f:
        reader = csv.reader(f, delimiter="\t")
        header = None
        gene_idx = -1
        status_idx = -1
        seq_idx = -1
        start_idx = -1
        end_idx = -1
        score_idx = -1

        for row in reader:
            if not row:
                continue
            first_val = row[0].strip()
            if not first_val:
                continue

            if header is None:
                header = [col.lstrip("#").strip().lower() for col in row]
                for idx, col in enumerate(header):
                    if col in ("gene", "gene_id", "busco_id"):
                        gene_idx = idx
                    elif col in ("status",):
                        status_idx = idx
                    elif col in ("sequence", "sequence_id", "seq_id", "seq"):
                        seq_idx = idx
                    elif col in ("gene start", "gene_start", "start"):
                        start_idx = idx
                    elif col in ("gene end", "gene_end", "end"):
                        end_idx = idx
                    elif col in ("score", "bitscore"):
                        score_idx = idx

                if gene_idx == -1:
                    gene_idx = 0
                if status_idx == -1:
                    status_idx = 1
                if seq_idx == -1:
                    seq_idx = 2
                if start_idx == -1:
                    start_idx = 3
                if end_idx == -1:
                    end_idx = 4
                if score_idx == -1:
                    score_idx = 6
                continue

            max_idx = max(gene_idx, status_idx, seq_idx, start_idx, end_idx, score_idx)
            if len(row) <= max_idx:
                continue

            status = row[status_idx].strip()
            if status.lower() == "missing":
                continue

            seq_id = row[seq_idx].strip()
            if not seq_id or seq_id in ("-", "n/a", "none"):
                continue

            gene_id = row[gene_idx].strip()
            start_s = row[start_idx].strip()
            end_s = row[end_idx].strip()
            score_s = row[score_idx].strip()

            start = int(float(start_s)) if start_s else 0
            end = int(float(end_s)) if end_s else 0
            score = float(score_s) if score_s else float(0)

            hit = CompleasmHit(
                gene_id=gene_id,
                status=status,
                sequence_id=seq_id,
                start=start,
                end=end,
                score=score,
            )
            hits.setdefault(seq_id, []).append(hit)

    return hits


def aggregate_node_annotations(
    tiara: dict[str, str | TiaraHit] | str | Path | None = None,
    tidk: dict[str, list[TidkWindow]] | str | Path | None = None,
    barrnap: dict[str, list[BarrnapHit]] | str | Path | None = None,
    compleasm: dict[str, list[CompleasmHit]] | str | Path | None = None,
    sequence_ids: Iterable[str] | None = None,
) -> dict[str, NodeAnnotations]:
    """Aggregate multi-layer annotations across unitig sequence IDs.

    Args:
        tiara: Tiara classification dictionary, file path, or None.
        tidk: TIDK windows dictionary, file path, or None.
        barrnap: Barrnap hits dictionary, file path, or None.
        compleasm: Compleasm hits dictionary, file path, or None.
        sequence_ids: Optional sequence IDs to initialize, ensuring presence in output.

    Returns:
        Mapping from sequence_id to aggregated NodeAnnotations record.
    """
    tiara_map: dict[str, str | TiaraHit] | None
    if isinstance(tiara, str | Path):
        tiara_map = load_tiara(tiara)
    else:
        tiara_map = tiara

    tidk_map: dict[str, list[TidkWindow]] | None
    if isinstance(tidk, str | Path):
        tidk_map = load_tidk(tidk)
    else:
        tidk_map = tidk

    barrnap_map: dict[str, list[BarrnapHit]] | None
    if isinstance(barrnap, str | Path):
        barrnap_map = load_barrnap(barrnap)
    else:
        barrnap_map = barrnap

    compleasm_map: dict[str, list[CompleasmHit]] | None
    if isinstance(compleasm, str | Path):
        compleasm_map = load_compleasm(compleasm)
    else:
        compleasm_map = compleasm

    all_ids: set[str] = set()
    if sequence_ids is not None:
        all_ids.update(sequence_ids)
    if tiara_map is not None:
        all_ids.update(tiara_map.keys())
    if tidk_map is not None:
        all_ids.update(tidk_map.keys())
    if barrnap_map is not None:
        all_ids.update(barrnap_map.keys())
    if compleasm_map is not None:
        all_ids.update(compleasm_map.keys())

    aggregated: dict[str, NodeAnnotations] = {}
    for seq_id in sorted(all_ids):
        t_val = tiara_map.get(seq_id) if tiara_map else None
        if isinstance(t_val, str):
            t_hit: TiaraHit | None = TiaraHit(sequence_id=seq_id, classification=t_val)
        elif isinstance(t_val, TiaraHit):
            t_hit = t_val
        else:
            t_hit = None

        w_hits = list(tidk_map.get(seq_id, [])) if tidk_map else []
        b_hits = list(barrnap_map.get(seq_id, [])) if barrnap_map else []
        c_hits = list(compleasm_map.get(seq_id, [])) if compleasm_map else []

        aggregated[seq_id] = NodeAnnotations(
            sequence_id=seq_id,
            tiara=t_hit,
            tidk_windows=w_hits,
            barrnap_hits=b_hits,
            compleasm_hits=c_hits,
        )

    return aggregated
