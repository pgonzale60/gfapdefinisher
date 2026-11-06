"""Out-of-core SQLite storage and indexer for Hifiasm A-line read placements.

Provides constant-memory streaming ingestion of GFA A-lines from uncompressed
or gzip-compressed files into an embedded SQLite database with post-ingestion
indexing and window/spanning query APIs.
"""

from __future__ import annotations

import gzip
import sqlite3
from pathlib import Path
from typing import Any, NamedTuple

DEFAULT_BATCH_SIZE: int = 100_000
MIN_A_LINE_TOKENS: int = 7
TAG_SPLIT_PARTS: int = 2


class ReadPlacement(NamedTuple):
    """Structured representation of a single GFA A-line read placement."""

    segment_name: str
    seg_pos: int
    seg_strand: str
    read_id: str
    read_pos: int
    read_len: int
    haplogroup: str | None
    raw_tags: str | None


class ReadPlacementStore:
    """Out-of-core SQLite store for GFA A-line read placements."""

    def __init__(self, db_path: str | Path | None = None) -> None:
        """Initialize the store backed by SQLite file or in-memory database."""
        if db_path is None or str(db_path) == ":memory:":
            self.db_path = ":memory:"
        else:
            self.db_path = str(db_path)
            p = Path(self.db_path)
            if p.parent and not p.parent.exists():
                p.parent.mkdir(parents=True, exist_ok=True)

        self._conn: sqlite3.Connection | None = sqlite3.connect(self.db_path)
        self._init_db()

    @property
    def conn(self) -> sqlite3.Connection:
        """Active SQLite connection handle."""
        if self._conn is None:
            raise RuntimeError("ReadPlacementStore connection is closed.")
        return self._conn

    def _init_db(self) -> None:
        """Initialize SQLite PRAGMAs and placement table schema."""
        self.conn.execute("PRAGMA synchronous = OFF")
        self.conn.execute("PRAGMA journal_mode = MEMORY")
        self.conn.execute(
            """
            CREATE TABLE IF NOT EXISTS read_placements (
                segment_name TEXT,
                seg_pos INTEGER,
                seg_strand TEXT,
                read_id TEXT,
                read_pos INTEGER,
                read_len INTEGER,
                haplogroup TEXT,
                raw_tags TEXT
            )
            """
        )
        self.conn.commit()

    def __enter__(self) -> ReadPlacementStore:
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_val: BaseException | None,
        exc_tb: Any,
    ) -> None:
        self.close()

    def close(self) -> None:
        """Commit outstanding transactions and close SQLite connection."""
        if self._conn is not None:
            self._conn.commit()
            self._conn.close()
            self._conn = None

    def ingest_gfa_a_lines(
        self,
        gfa_path: str | Path,
        batch_size: int = DEFAULT_BATCH_SIZE,
    ) -> int:
        """Stream and ingest GFA A-line read placements without heap exhaustion.

        Streams text or gzip-compressed GFA files line-by-line, parsing A-lines:
            A <seg_name> <seg_pos> <seg_strand> <read_id> <read_pos> <read_len> [tags...]
        and inserts in batches using executemany with PRAGMA synchronous = OFF and
        PRAGMA journal_mode = MEMORY. Indices on (segment_name, seg_pos) and read_id
        are created after the bulk insert completes.

        Args:
            gfa_path: Path to plain text or gzip (.gz) GFA file.
            batch_size: Batch size for executemany inserts.

        Returns:
            Count of ingested A-line records.
        """
        path = Path(gfa_path)
        if not path.is_file():
            raise FileNotFoundError(f"GFA file not found: {path}")

        is_gz = False
        if path.suffix == ".gz":
            is_gz = True
        else:
            with open(path, "rb") as probe:
                magic = probe.read(2)
                if magic == b"\x1f\x8b":
                    is_gz = True

        opener = gzip.open if is_gz else open

        count = 0
        batch: list[tuple[str, int, str, str, int, int, str | None, str | None]] = []
        cursor = self.conn.cursor()

        insert_sql = (
            "INSERT INTO read_placements ("
            "segment_name, seg_pos, seg_strand, read_id, "
            "read_pos, read_len, haplogroup, raw_tags"
            ") VALUES (?, ?, ?, ?, ?, ?, ?, ?)"
        )

        with opener(path, "rt", encoding="utf-8", errors="replace") as stream:
            for line in stream:
                if not line.startswith("A\t"):
                    continue

                tokens = line.rstrip("\r\n").split("\t")
                if len(tokens) < MIN_A_LINE_TOKENS:
                    continue

                (
                    _rec_type,
                    seg_name,
                    seg_pos_str,
                    seg_strand,
                    read_id,
                    read_pos_str,
                    read_len_str,
                    *tags,
                ) = tokens

                haplogroup: str | None = None
                raw_tags: str | None = None

                if tags:
                    raw_tags = "\t".join(tags)
                    for tag in tags:
                        if tag.startswith("HG:"):
                            parts = tag.split(":", TAG_SPLIT_PARTS)
                            haplogroup = parts[-1]
                            break
                    if haplogroup is None and ":" not in tags[0]:
                        haplogroup = tags[0]

                batch.append(
                    (
                        seg_name,
                        int(seg_pos_str),
                        seg_strand,
                        read_id,
                        int(read_pos_str),
                        int(read_len_str),
                        haplogroup,
                        raw_tags,
                    )
                )

                if len(batch) >= batch_size:
                    cursor.executemany(insert_sql, batch)
                    count += len(batch)
                    batch.clear()

        if batch:
            cursor.executemany(insert_sql, batch)
            count += len(batch)
            batch.clear()

        self.conn.commit()

        # Build indices *after* bulk insert
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_read_placements_seg ON read_placements (segment_name, seg_pos)")
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_read_placements_read ON read_placements (read_id)")
        self.conn.commit()

        return count

    def get_reads_for_segment(
        self,
        segment_name: str,
        start: int = 0,
        end: int = -1,
    ) -> list[tuple]:
        """Retrieve read placements for a given segment within an optional window.

        Args:
            segment_name: Segment identifier.
            start: 0-based coordinate start bound (inclusive, default 0).
            end: Upper coordinate bound (inclusive, or -1 for unbound).

        Returns:
            List of read placement tuples ordered by seg_pos.
        """
        cursor = self.conn.cursor()
        if end < 0:
            if start <= 0:
                query = (
                    "SELECT segment_name, seg_pos, seg_strand, read_id, "
                    "read_pos, read_len, haplogroup, raw_tags "
                    "FROM read_placements "
                    "WHERE segment_name = ? "
                    "ORDER BY seg_pos"
                )
                cursor.execute(query, (segment_name,))
            else:
                query = (
                    "SELECT segment_name, seg_pos, seg_strand, read_id, "
                    "read_pos, read_len, haplogroup, raw_tags "
                    "FROM read_placements "
                    "WHERE segment_name = ? AND seg_pos >= ? "
                    "ORDER BY seg_pos"
                )
                cursor.execute(query, (segment_name, start))
        else:
            query = (
                "SELECT segment_name, seg_pos, seg_strand, read_id, "
                "read_pos, read_len, haplogroup, raw_tags "
                "FROM read_placements "
                "WHERE segment_name = ? AND seg_pos >= ? AND seg_pos <= ? "
                "ORDER BY seg_pos"
            )
            cursor.execute(query, (segment_name, start, end))

        return cursor.fetchall()

    def get_spanning_reads(self, seg1: str, seg2: str) -> list[str]:
        """Retrieve distinct read IDs that map to both seg1 and seg2.

        Args:
            seg1: First segment name.
            seg2: Second segment name.

        Returns:
            List of distinct read IDs spanning both segments, sorted lexicographically.
        """
        cursor = self.conn.cursor()
        query = (
            "SELECT read_id FROM read_placements WHERE segment_name = ? "
            "INTERSECT "
            "SELECT read_id FROM read_placements WHERE segment_name = ? "
            "ORDER BY read_id"
        )
        cursor.execute(query, (seg1, seg2))
        return [row[0] for row in cursor.fetchall()]

    def get_read_count(self) -> int:
        """Return total number of read placement records in the store."""
        cursor = self.conn.cursor()
        cursor.execute("SELECT COUNT(*) FROM read_placements")
        row = cursor.fetchone()
        return int(row[0]) if row is not None else 0

    def get_distinct_read_count(self) -> int:
        """Return total number of distinct read IDs in the store."""
        cursor = self.conn.cursor()
        cursor.execute("SELECT COUNT(DISTINCT read_id) FROM read_placements")
        row = cursor.fetchone()
        return int(row[0]) if row is not None else 0

    def __len__(self) -> int:
        """Return total placement count."""
        return self.get_read_count()
