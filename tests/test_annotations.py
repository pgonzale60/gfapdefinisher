"""Unit and integration tests for precomputed annotation loaders."""

import gzip
from pathlib import Path

import pytest

from gfapdefinisher.annotations import (
    BarrnapHit,
    CompleasmHit,
    NodeAnnotations,
    TiaraHit,
    TidkWindow,
    aggregate_node_annotations,
    load_barrnap,
    load_compleasm,
    load_tiara,
    load_tidk,
)

BENCHMARK_DIR = Path("data/benchmark/DF5120")


def test_data_structures():
    """Verify dataclass slots, immutability, and attribute semantics."""
    tiara = TiaraHit(sequence_id="utg1", classification="eukarya")
    assert tiara.sequence_id == "utg1"
    assert tiara.classification == "eukarya"
    with pytest.raises(AttributeError):
        tiara.classification = "bacteria"  # type: ignore[misc]

    tidk = TidkWindow(sequence_id="utg1", window_start=1000, forward_count=5, reverse_count=0)
    assert tidk.sequence_id == "utg1"
    assert tidk.window_start == 1000
    assert tidk.forward_count == 5
    assert tidk.reverse_count == 0

    barrnap = BarrnapHit(
        sequence_id="utg1",
        start=100,
        end=200,
        strand="+",
        feature_name="5S_rRNA",
        evalue=1e-10,
    )
    assert barrnap.feature_name == "5S_rRNA"
    assert barrnap.evalue == 1e-10

    compleasm = CompleasmHit(
        gene_id="12345at119089",
        status="Single",
        sequence_id="utg1",
        start=5000,
        end=9000,
        score=1500.5,
    )
    assert compleasm.gene_id == "12345at119089"
    assert compleasm.status == "Single"
    assert compleasm.score == 1500.5


def test_node_annotations_container():
    """Verify NodeAnnotations container slots, properties, and helper methods."""
    node = NodeAnnotations(sequence_id="utg1")
    assert node.sequence_id == "utg1"
    assert node.tiara is None
    assert node.tiara_classification is None
    assert node.tiara_hit is None
    assert not node.has_telomere()
    assert not node.has_rdna()
    assert not node.is_eukaryotic()
    assert not node.is_contaminant()

    # String initialization for tiara converts to TiaraHit
    node_str_tiara = NodeAnnotations(sequence_id="utg2", tiara="eukarya")  # type: ignore[arg-type]
    assert node_str_tiara.tiara == TiaraHit(sequence_id="utg2", classification="eukarya")
    assert node_str_tiara.tiara_classification == "eukarya"
    assert node_str_tiara.is_eukaryotic()
    assert not node_str_tiara.is_contaminant()

    # Bacteria classification
    node_bact = NodeAnnotations(
        sequence_id="utg3",
        tiara=TiaraHit(sequence_id="utg3", classification="bacteria"),
    )
    assert node_bact.is_contaminant()
    assert not node_bact.is_eukaryotic()

    # TIDK, Barrnap, Compleasm aggregation
    tidk_win = TidkWindow(sequence_id="utg1", window_start=1000, forward_count=10, reverse_count=0)
    barrnap_hit = BarrnapHit(
        sequence_id="utg1",
        start=500,
        end=600,
        strand="-",
        feature_name="18S_rRNA",
        evalue=2e-20,
    )
    compleasm_hit = CompleasmHit(
        gene_id="g1",
        status="Single",
        sequence_id="utg1",
        start=1000,
        end=2000,
        score=500.0,
    )

    node.tidk_windows.append(tidk_win)
    node.barrnap_hits.append(barrnap_hit)
    node.compleasm_hits.append(compleasm_hit)

    assert node.has_telomere(threshold=5)
    assert not node.has_telomere(threshold=20)
    assert node.has_rdna()
    assert node.tidk == [tidk_win]
    assert node.barrnap == [barrnap_hit]
    assert node.compleasm == [compleasm_hit]


def test_load_tiara_mock(tmp_path: Path):
    """Verify load_tiara parser with mock TSV and gzipped TSV."""
    content = (
        "# Tiara mock output\n"
        "sequence_id\tclass_fst_stage\tclass_snd_stage\torg\tbac\tarc\teuk\tunk1\tpla\tunk2\tmit\n"
        "utg001\teukarya\tn/a\t0.01\t0.02\t0.01\t0.96\t0.00\tn/a\tn/a\tn/a\n"
        "utg002\tbacteria\tn/a\t0.00\t1.00\t0.00\t0.00\t0.00\tn/a\tn/a\tn/a\n"
        "utg003\torganelle\tplastid\t0.90\t0.05\t0.00\t0.05\t0.00\t0.90\t0.00\t0.10\n"
    )
    tsv_path = tmp_path / "tiara.txt"
    tsv_path.write_text(content, encoding="utf-8")

    result = load_tiara(tsv_path)
    assert result == {
        "utg001": "eukarya",
        "utg002": "bacteria",
        "utg003": "organelle",
    }

    # Gzip test
    gz_path = tmp_path / "tiara.txt.gz"
    with gzip.open(gz_path, "wt", encoding="utf-8") as f:
        f.write(content)

    gz_result = load_tiara(gz_path)
    assert gz_result == result


def test_load_tidk_mock(tmp_path: Path):
    """Verify load_tidk parser with mock TSV."""
    content = (
        "id\twindow\tforward_repeat_number\treverse_repeat_number\ttelomeric_repeat\n"
        "utg001\t1000\t0\t0\tTTAGGC\n"
        "utg001\t2000\t12\t1\tTTAGGC\n"
        "utg002\t1000\t3\t15\tTTAGGC\n"
    )
    tsv_path = tmp_path / "tidk.tsv"
    tsv_path.write_text(content, encoding="utf-8")

    result = load_tidk(tsv_path)
    assert len(result) == 2
    assert "utg001" in result
    assert "utg002" in result

    assert len(result["utg001"]) == 2
    assert result["utg001"][0] == TidkWindow("utg001", 1000, 0, 0)
    assert result["utg001"][1] == TidkWindow("utg001", 2000, 12, 1)

    assert len(result["utg002"]) == 1
    assert result["utg002"][0] == TidkWindow("utg002", 1000, 3, 15)


def test_load_barrnap_mock(tmp_path: Path):
    """Verify load_barrnap parser with mock GFF3 annotations."""
    content = (
        "##gff-version 3\n"
        "# Comment line\n"
        "utg001\tbarrnap:0.9\trRNA\t100\t200\t3.5e-16\t+\t.\tName=5S_rRNA;product=5S ribosomal RNA\n"
        "utg001\tbarrnap:0.9\trRNA\t500\t1500\t1.2e-40\t-\t.\tName=18S_rRNA\n"
        "utg002\tbarrnap:0.9\trRNA\t300\t400\t.\t+\t.\tID=custom_rrna;product=28S ribosomal RNA\n"
    )
    gff_path = tmp_path / "barrnap.gff"
    gff_path.write_text(content, encoding="utf-8")

    result = load_barrnap(gff_path)
    assert len(result) == 2
    assert len(result["utg001"]) == 2
    assert result["utg001"][0] == BarrnapHit(
        sequence_id="utg001",
        start=100,
        end=200,
        strand="+",
        feature_name="5S_rRNA",
        evalue=3.5e-16,
    )
    assert result["utg001"][1] == BarrnapHit(
        sequence_id="utg001",
        start=500,
        end=1500,
        strand="-",
        feature_name="18S_rRNA",
        evalue=1.2e-40,
    )

    assert len(result["utg002"]) == 1
    # Fallback to ID or product
    assert result["utg002"][0].feature_name == "custom_rrna"
    assert result["utg002"][0].evalue == 0.0


def test_load_compleasm_mock(tmp_path: Path):
    """Verify load_compleasm parser with mock USCO table."""
    content = (
        "Gene\tStatus\tSequence\tGene Start\tGene End\tStrand\tScore\tLength\tIdentity\tFraction\tFrameshift events\tBest gene\tCodons\n"
        "gene1\tSingle\tutg001\t1000\t2000\t+\t500.5\t300\t0.9\t0.8\t0\tg1_best\t...\n"
        "gene2\tDuplicated\tutg001\t3000\t4500\t-\t1200\t400\t0.95\t0.9\t1\tg2_best\t...\n"
        "gene3\tMissing\n"
        "gene4\tFragmented\tutg002\t500\t800\t+\t250\t100\t0.5\t0.4\t0\tg4_best\t...\n"
    )
    tsv_path = tmp_path / "compleasm.tsv"
    tsv_path.write_text(content, encoding="utf-8")

    result = load_compleasm(tsv_path)
    assert len(result) == 2
    assert "utg001" in result
    assert "utg002" in result

    assert len(result["utg001"]) == 2
    assert result["utg001"][0] == CompleasmHit(
        gene_id="gene1",
        status="Single",
        sequence_id="utg001",
        start=1000,
        end=2000,
        score=500.5,
    )
    assert result["utg001"][1] == CompleasmHit(
        gene_id="gene2",
        status="Duplicated",
        sequence_id="utg001",
        start=3000,
        end=4500,
        score=1200.0,
    )

    assert len(result["utg002"]) == 1
    assert result["utg002"][0] == CompleasmHit(
        gene_id="gene4",
        status="Fragmented",
        sequence_id="utg002",
        start=500,
        end=800,
        score=250.0,
    )


def test_aggregate_node_annotations_mock():
    """Verify merging and union of multi-layer annotation dictionaries."""
    tiara = {"utg1": "eukarya", "utg2": "bacteria"}
    tidk = {
        "utg1": [TidkWindow("utg1", 1000, 15, 0)],
        "utg3": [TidkWindow("utg3", 500, 2, 2)],
    }
    barrnap = {
        "utg1": [BarrnapHit("utg1", 10, 50, "+", "5S_rRNA", 1e-5)],
    }
    compleasm = {
        "utg1": [CompleasmHit("g1", "Single", "utg1", 100, 200, 300.0)],
        "utg4": [CompleasmHit("g2", "Single", "utg4", 500, 600, 400.0)],
    }

    # Sequence IDs explicitly specifies utg5 which has no hits in any layer
    agg = aggregate_node_annotations(
        tiara=tiara,
        tidk=tidk,
        barrnap=barrnap,
        compleasm=compleasm,
        sequence_ids=["utg1", "utg5"],
    )

    assert set(agg.keys()) == {"utg1", "utg2", "utg3", "utg4", "utg5"}

    # utg1 has hits in all 4 layers
    u1 = agg["utg1"]
    assert u1.is_eukaryotic()
    assert u1.has_telomere()
    assert u1.has_rdna()
    assert len(u1.compleasm_hits) == 1

    # utg2 only has tiara
    u2 = agg["utg2"]
    assert u2.is_contaminant()
    assert not u2.has_telomere()
    assert not u2.has_rdna()
    assert len(u2.compleasm_hits) == 0

    # utg5 has empty annotations
    u5 = agg["utg5"]
    assert u5.tiara is None
    assert len(u5.tidk_windows) == 0
    assert len(u5.barrnap_hits) == 0
    assert len(u5.compleasm_hits) == 0


@pytest.mark.benchmark
def test_benchmark_df5120_staged_integration():
    """Integration test parsing staged benchmark files if present in workspace."""
    if not BENCHMARK_DIR.exists():
        pytest.skip(f"Benchmark directory {BENCHMARK_DIR} not present.")

    tiara_path = BENCHMARK_DIR / "tiara_unitigs.txt"
    tidk_path = BENCHMARK_DIR / "tidk_telomeres_telomeric_repeat_windows.tsv"
    barrnap_path = BENCHMARK_DIR / "barrnap_rdna.gff"
    compleasm_path = BENCHMARK_DIR / "compleasm_unitigs_full_table.tsv"

    assert tiara_path.exists()
    assert tidk_path.exists()
    assert barrnap_path.exists()
    assert compleasm_path.exists()

    tiara = load_tiara(tiara_path)
    assert len(tiara) > 200
    assert "utg000001l" in tiara
    assert tiara["utg000001l"] == "eukarya"

    tidk = load_tidk(tidk_path)
    assert len(tidk) > 200
    assert "utg000001l" in tidk
    assert len(tidk["utg000001l"]) > 5000

    barrnap = load_barrnap(barrnap_path)
    assert len(barrnap) > 30
    assert "utg000001l" in barrnap
    assert any("5S_rRNA" in h.feature_name for h in barrnap["utg000001l"])

    compleasm = load_compleasm(compleasm_path)
    assert len(compleasm) > 15
    assert "utg000019l" in compleasm
    assert any(h.status == "Single" for h in compleasm["utg000019l"])

    # Aggregate by paths directly
    agg = aggregate_node_annotations(
        tiara=tiara_path,
        tidk=tidk_path,
        barrnap=barrnap_path,
        compleasm=compleasm_path,
    )
    assert len(agg) >= 222
    assert "utg000001l" in agg
    assert agg["utg000001l"].is_eukaryotic()
    assert agg["utg000001l"].has_telomere()
    assert agg["utg000001l"].has_rdna()
