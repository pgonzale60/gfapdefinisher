"""Unit and integration tests for ReadPlacementStore out-of-core SQLite storage."""

from __future__ import annotations

import gzip
import sqlite3
from pathlib import Path

import pytest

from gfapdefinisher.store import ReadPlacement, ReadPlacementStore


def test_in_memory_initialization_and_context_manager():
    """Verify in-memory SQLite store initialization and context manager protocol."""
    with ReadPlacementStore() as store:
        assert store.get_read_count() == 0
        assert store.get_distinct_read_count() == 0
        assert len(store) == 0
        assert isinstance(store.conn, sqlite3.Connection)

    # After exit, connection is closed
    with pytest.raises(RuntimeError, match="closed"):
        _ = store.conn


def test_file_backed_store_lifecycle(tmp_path: Path):
    """Verify file-backed SQLite store creation, persistence, and teardown."""
    db_file = tmp_path / "test_store.db"
    assert not db_file.exists()

    with ReadPlacementStore(db_file) as store:
        assert db_file.exists()
        assert store.get_read_count() == 0

    assert db_file.exists()
    # Re-open existing db
    with ReadPlacementStore(db_file) as store:
        assert store.get_read_count() == 0


def test_synthetic_gfa_ingestion_and_queries(tmp_path: Path):
    """Test streaming ingestion, count, indexing, and window/spanning query APIs."""
    gfa_content = "\n".join(
        [
            "H\tVN:Z:1.0",
            "S\tutg0001\t1000\t*\tLN:i:1000",
            "S\tutg0002\t2000\t*\tLN:i:2000",
            "S\tutg0003\t1500\t*\tLN:i:1500",
            "# Comment line should be ignored",
            "",
            # Read 1: on utg0001 (pos 50) and utg0002 (pos 200) -> spans utg0001 & utg0002
            "A\tutg0001\t50\t+\tread_alpha\t0\t15000\tHG:A:p\tNM:i:0",
            "A\tutg0002\t200\t-\tread_alpha\t1000\t15000\tHG:A:p\tNM:i:1",
            # Read 2: on utg0001 (pos 150) and utg0002 (pos 350) -> spans utg0001 & utg0002
            "A\tutg0001\t150\t+\tread_beta\t0\t18000\tHG:A:m",
            "A\tutg0002\t350\t+\tread_beta\t2000\t18000\tHG:A:m",
            # Read 3: only on utg0001 (pos 500)
            "A\tutg0001\t500\t+\tread_gamma\t100\t12000\tp",
            # Read 4: only on utg0003
            "A\tutg0003\t100\t+\tread_delta\t0\t10000",
            # Malformed short A-line (should be skipped safely)
            "A\tshort_line",
            # Dovetail link line
            "L\tutg0001\t+\tutg0002\t+\t500M",
        ]
    )

    gfa_path = tmp_path / "test.gfa"
    gfa_path.write_text(gfa_content, encoding="utf-8")

    db_path = tmp_path / "placements.db"
    with ReadPlacementStore(db_path) as store:
        # Ingest with tiny batch size to test batch boundaries and post-insert indexing
        ingested = store.ingest_gfa_a_lines(gfa_path, batch_size=2)
        assert ingested == 6
        assert store.get_read_count() == 6
        assert store.get_distinct_read_count() == 4
        assert len(store) == 6

        # Query segment utg0001 without window
        reads_u1 = store.get_reads_for_segment("utg0001")
        assert len(reads_u1) == 3
        # Ordered by seg_pos: 50, 150, 500
        assert reads_u1[0][0] == "utg0001"
        assert reads_u1[0][1] == 50
        assert reads_u1[0][3] == "read_alpha"
        assert reads_u1[0][6] == "p"
        assert reads_u1[0][7] == "HG:A:p\tNM:i:0"

        assert reads_u1[1][1] == 150
        assert reads_u1[1][3] == "read_beta"
        assert reads_u1[1][6] == "m"

        assert reads_u1[2][1] == 500
        assert reads_u1[2][3] == "read_gamma"
        assert reads_u1[2][6] == "p"

        # Query segment utg0001 with window bounds
        window_reads = store.get_reads_for_segment("utg0001", start=100, end=400)
        assert len(window_reads) == 1
        assert window_reads[0][3] == "read_beta"

        # Query segment utg0001 with start bound only
        start_reads = store.get_reads_for_segment("utg0001", start=100, end=-1)
        assert len(start_reads) == 2
        assert [r[3] for r in start_reads] == ["read_beta", "read_gamma"]

        # Non-existent segment
        empty_reads = store.get_reads_for_segment("utg_missing")
        assert empty_reads == []

        # Spanning reads between utg0001 and utg0002
        spanning = store.get_spanning_reads("utg0001", "utg0002")
        assert spanning == ["read_alpha", "read_beta"]

        # Spanning reads between utg0001 and utg0003 (no common reads)
        spanning_none = store.get_spanning_reads("utg0001", "utg0003")
        assert spanning_none == []

        # Verify index presence in sqlite master
        indices = store.conn.execute(
            "SELECT name FROM sqlite_master WHERE type='index' AND name LIKE 'idx_read_placements%'"
        ).fetchall()
        idx_names = {row[0] for row in indices}
        assert "idx_read_placements_seg" in idx_names
        assert "idx_read_placements_read" in idx_names

        # Check ReadPlacement helper
        first_parsed = ReadPlacement(*reads_u1[0])
        assert first_parsed.segment_name == "utg0001"
        assert first_parsed.seg_pos == 50
        assert first_parsed.read_id == "read_alpha"
        assert first_parsed.haplogroup == "p"


def test_gzip_gfa_ingestion(tmp_path: Path):
    """Test ingestion from gzip-compressed GFA files."""
    gfa_lines = [
        "A\tutg1\t10\t+\tr1\t0\t1000\tHG:A:p",
        "A\tutg2\t20\t-\tr2\t0\t1000\tHG:A:m",
    ]
    gz_path = tmp_path / "sample.gfa.gz"
    with gzip.open(gz_path, "wt", encoding="utf-8") as gz_out:
        gz_out.write("\n".join(gfa_lines) + "\n")

    with ReadPlacementStore() as store:
        count = store.ingest_gfa_a_lines(gz_path)
        assert count == 2
        assert store.get_read_count() == 2
        reads = store.get_reads_for_segment("utg1")
        assert len(reads) == 1
        assert reads[0][3] == "r1"


def test_file_not_found(tmp_path: Path):
    """Test that missing GFA file raises FileNotFoundError."""
    missing = tmp_path / "does_not_exist.gfa"
    with ReadPlacementStore() as store:
        with pytest.raises(FileNotFoundError, match="not found"):
            store.ingest_gfa_a_lines(missing)


@pytest.mark.benchmark
def test_benchmark_data_streaming_if_present():
    """Integration test verifying constant-memory streaming against benchmark data if present."""
    bench_file = Path("data") / "benchmark" / "DF5120" / "asm.p_utg.gfa.gz"
    if not bench_file.exists():
        pytest.skip(f"Benchmark file {bench_file} not present in environment")

    with ReadPlacementStore() as store:
        ingested = store.ingest_gfa_a_lines(bench_file, batch_size=10000)
        assert ingested > 0
        assert store.get_read_count() == ingested


def test_empty_gfa_file(tmp_path: Path):
    """Test ingestion from an empty GFA file."""
    empty_file = tmp_path / "empty.gfa"
    empty_file.write_text("", encoding="utf-8")

    with ReadPlacementStore() as store:
        count = store.ingest_gfa_a_lines(empty_file)
        assert count == 0
        assert store.get_read_count() == 0
        assert store.get_spanning_reads("s1", "s2") == []
        assert store.get_reads_for_segment("s1") == []


def test_large_batch_streaming(tmp_path: Path):
    """Test streaming behavior across multiple batch boundaries."""
    total_lines = 250
    batch_size = 50
    lines = [f"A\tseg_{i % 5}\t{i * 10}\t+\tread_{i}\t0\t15000\tHG:A:p" for i in range(total_lines)]
    gfa_file = tmp_path / "large.gfa"
    gfa_file.write_text("\n".join(lines) + "\n", encoding="utf-8")

    with ReadPlacementStore() as store:
        ingested = store.ingest_gfa_a_lines(gfa_file, batch_size=batch_size)
        assert ingested == total_lines
        assert store.get_read_count() == total_lines
        assert store.get_distinct_read_count() == total_lines

        # Seg 0 should have 50 reads
        seg0_reads = store.get_reads_for_segment("seg_0")
        assert len(seg0_reads) == 50


def test_duplicate_reads_and_spanning_deduplication(tmp_path: Path):
    """Test that reads appearing multiple times on segments are deduplicated in spanning query."""
    lines = [
        "A\tsegA\t100\t+\trShared\t0\t10000",
        "A\tsegA\t200\t+\trShared\t500\t10000",  # same read on segA twice
        "A\tsegB\t300\t-\trShared\t1000\t10000",
        "A\tsegB\t400\t-\trShared\t1500\t10000",  # same read on segB twice
        "A\tsegA\t500\t+\trPrivateA\t0\t10000",
        "A\tsegB\t600\t+\trPrivateB\t0\t10000",
    ]
    gfa_file = tmp_path / "multi.gfa"
    gfa_file.write_text("\n".join(lines) + "\n", encoding="utf-8")

    with ReadPlacementStore() as store:
        store.ingest_gfa_a_lines(gfa_file)
        assert store.get_read_count() == 6
        assert store.get_distinct_read_count() == 3

        spanning = store.get_spanning_reads("segA", "segB")
        assert spanning == ["rShared"]
