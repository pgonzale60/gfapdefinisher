"""Reproduce reviewed chromosome paths and read bridges from a curation recipe."""

from __future__ import annotations

import hashlib
import json
import os
import re
from collections import defaultdict
from io import DEFAULT_BUFFER_SIZE
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Any

from gfapdefinisher.emitter import emit_fasta
from gfapdefinisher.graph import load_gfa
from gfapdefinisher.policy import DEFAULT_POLICY
from gfapdefinisher.resolver import load_reads_sequences, reverse_complement
from gfapdefinisher.telomere import scan_terminal_tract

_GAP = re.compile("N+")
_GRAPH_SOURCES = {"primary", "raw"}
_READ_SOURCES = {"hifi", "ont"}
_GAP_TYPES = {"rdna", "satellite_port"}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(DEFAULT_BUFFER_SIZE), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _oriented(sequences: dict[str, str], reference: dict[str, Any]) -> str:
    name = reference["record"]
    if name not in sequences or not sequences[name]:
        raise ValueError(f"Missing sequence for {name}")
    seq = sequences[name].upper()
    strand = reference.get("strand", "+")
    if strand == "+":
        return seq
    if strand == "-":
        return reverse_complement(seq)
    raise ValueError(f"{name}: invalid strand {strand!r}")


def _bounds(reference: dict[str, Any], length: int) -> tuple[int, int]:
    start, end = reference.get("start", 0), reference.get("end", length)
    if not 0 <= start < end <= length:
        raise ValueError(f"{reference['record']}: slice [{start}, {end}) outside {length} bp sequence")
    return start, end


def _assemble_path(parts: list[dict[str, Any]], sources: dict[str, dict[str, str]], gap_bp: int, anchor_bp: int) -> str:
    oriented: list[tuple[str, int, int]] = []
    for part in parts:
        if "gap" in part:
            if part["gap"] not in _GAP_TYPES:
                raise ValueError(f"Unknown gap type {part['gap']!r}")
            oriented.append(("N" * gap_bp, 0, gap_bp))
        else:
            seq = _oriented(sources[part["source"]], part)
            start, end = _bounds(part, len(seq))
            oriented.append((seq, start, end))
    pieces = [seq[start:end] for seq, start, end in oriented]
    for index, part in enumerate(parts):
        seq, start, end = oriented[index]
        source = part.get("source")
        if source in _READ_SOURCES:
            if not 0 < index < len(parts) - 1 or not anchor_bp <= start < end <= len(seq) - anchor_bp:
                raise ValueError(f"{part['record']}: read bridge lacks two flanks")
            left = pieces[index - 1][-anchor_bp:]
            right = pieces[index + 1][:anchor_bp]
            if left != seq[start - anchor_bp : start] or right != seq[end : end + anchor_bp]:
                raise ValueError(f"{part['record']}: read bridge does not match both splice anchors")
            if set(left + pieces[index] + right) - set("ACGT"):
                raise ValueError(f"{part['record']}: read bridge or its anchors contain ambiguous bases")
        elif (
            index
            and "gap" not in part
            and "gap" not in parts[index - 1]
            and parts[index - 1].get("source") not in _READ_SOURCES
        ):
            overlap = part["overlap_bp"]
            if not 0 < overlap <= min(start, len(pieces[index - 1])):
                raise ValueError(f"{part['record']}: invalid upstream overlap length")
            if pieces[index - 1][-overlap:] != seq[start - overlap : start]:
                raise ValueError(f"{part['record']}: upstream overlap does not match")
    return "".join(pieces)


def curate(
    recipe_path: Path,
    source_paths: dict[str, Path | None],
    output: Path,
) -> str:
    """Apply a reviewed recipe and return the SHA-256 of the validated FASTA.

    Recipe coordinates are zero-based and half-open in the oriented sequence.
    An existing output is reusable if its bytes match the validated result.
    """
    recipe = json.loads(recipe_path.read_text(encoding="utf-8"))
    if recipe["version"] != 1:
        raise ValueError("Unsupported curation recipe version")

    gap_bp, anchor_bp, cap_bp, line_width = (
        recipe["gap_bp"],
        recipe["bridge_anchor_bp"],
        recipe["cap_match_bp"],
        recipe["line_width"],
    )
    if min(gap_bp, anchor_bp, cap_bp, line_width) <= 0:
        raise ValueError("Gap, anchor, cap and FASTA line lengths must be positive")
    chromosomes = recipe["chromosomes"]
    names = [chrom["name"] for chrom in chromosomes]
    if (
        not names
        or len(set(names)) != len(names)
        or any(len(name.split()) != 1 or name != name.strip() for name in names)
    ):
        raise ValueError("The recipe must declare distinct, nonempty FASTA record IDs")
    inputs = {path.resolve() for path in source_paths.values() if path is not None}
    if output.resolve() in inputs | {recipe_path.resolve()}:
        raise ValueError("The output must differ from the inputs")

    required: dict[str, set[str]] = defaultdict(set)
    for chrom in chromosomes:
        ends = (chrom["cap_5p"], chrom["cap_3p"])
        if any(cap["source"] not in _GRAPH_SOURCES for cap in ends):
            raise ValueError("Assigned telomere caps must come from a graph")
        for item in (*chrom["parts"], *ends):
            if "gap" not in item:
                required[item["source"]].add(item["record"])
    sources: dict[str, dict[str, str]] = {}
    for label, records in required.items():
        path = source_paths.get(label)
        if path is None or not path.is_file():
            raise ValueError(f"Missing input for recipe source {label!r}")
        if label in _GRAPH_SOURCES:
            graph = load_gfa(path, load_sequences=True)
            sources[label] = {name: graph[name].sequence or "" for name in records if name in graph}
            del graph
        else:
            sources[label] = load_reads_sequences(path, records)

    result: dict[str, str] = {}
    for chrom in chromosomes:
        name = chrom["name"]
        sequence = _assemble_path(chrom["parts"], sources, gap_bp, anchor_bp)
        if set(sequence) - set("ACGTN"):
            raise ValueError(f"{name}: unexpected characters in curated sequence")
        if len(_GAP.findall(sequence)) != chrom["gap_count"]:
            raise ValueError(f"{name}: unexpected number of remaining gaps")
        sequence = _GAP.sub("N" * gap_bp, sequence)
        cap_5p = _oriented(sources[chrom["cap_5p"]["source"]], chrom["cap_5p"])
        cap_3p = _oriented(sources[chrom["cap_3p"]["source"]], chrom["cap_3p"])
        if (
            min(len(sequence), len(cap_5p), len(cap_3p)) < cap_bp
            or sequence[:cap_bp] != cap_5p[:cap_bp]
            or sequence[-cap_bp:] != cap_3p[-cap_bp:]
        ):
            raise ValueError(f"{name}: sequence ends do not match the assigned caps")
        telomere = DEFAULT_POLICY.telomere
        if not (
            scan_terminal_tract(sequence, telomere.motif_reverse.value).is_terminal_cap
            and scan_terminal_tract(sequence, telomere.motif_forward.value, is_suffix=True).is_terminal_cap
        ):
            raise ValueError(f"{name}: assigned caps lack terminal telomere tracts")
        result[name] = sequence

    output.parent.mkdir(parents=True, exist_ok=True)
    with NamedTemporaryFile(mode="w", encoding="utf-8", newline="\n", dir=output.parent) as temporary:
        emit_fasta(result, temporary, line_width=line_width)
        temporary.flush()
        digest = _sha256(Path(temporary.name))
        expected = recipe.get("expected_sha256")
        if expected is not None and digest != expected:
            raise ValueError(f"Curated FASTA SHA-256 is {digest}, expected {expected}")
        try:
            os.link(temporary.name, output)
        except FileExistsError:
            if _sha256(output) != digest:
                raise
    return digest
