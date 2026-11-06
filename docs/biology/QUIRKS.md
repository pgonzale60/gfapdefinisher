# Biological Quirk and End-State Taxonomy

This document establishes the ratified 3-layer biological feature taxonomy for `gfapdefinisher` and downstream consumers (e.g. `PDE_origin`).
Rather than silently flattening complex assembly graph features or collapsing variable repeats into generic `N` gaps, `gfapdefinisher` emits a typed biological ledger (`quirks.tsv`) and AGP uncertainty metadata alongside finished T2T pseudomolecules.

Ratified in [GitHub Issue #26](https://github.com/pgonzale60/gfapdefinisher/issues/26).

---

## 1. Motivation and Architectural Overview

Real assembly graphs deviate from idealized linear topologies in biologically significant ways:
- **Variable-length subterminal satellites**: Tandem satellite arrays whose repeat unit copy numbers differ among individual ONT reads span between core somatic arm unitigs and distal germline telomeric caps.
- **Asymmetric lineage preservation**: Programmed DNA Elimination (PDE) creates somatic diminution forks where somatic stubs may drop out while germline caps persist, or vice versa.
- **Lineage exemption**: Certain chromosome ends (such as the rDNA array on nematode Chr II) are exempt from diminution, sharing a single terminal telomere cap across somatic and germline lineages.

To preserve these biological truths without compromising T2T assembly completeness or producing unannotated `N` gaps (BIO-001, BIO-004), `gfapdefinisher` implements a 3-layer hierarchy:

```
┌─────────────────────────────────────────────────────────────┐
│ Layer 1: Chromosome End State (Scope: end)                  │
│ Evaluates both lineages for 1 chromosome arm terminus       │
└──────────────────────────────┬──────────────────────────────┘
                               │
┌──────────────────────────────▼──────────────────────────────┐
│ Layer 2: Diminution Junction (Scope: junction)              │
│ Classifies each somatic diminution fork or break point      │
└──────────────────────────────┬──────────────────────────────┘
                               │
┌──────────────────────────────▼──────────────────────────────┐
│ Layer 3: Interval Biological Quirks (Scope: interval)        │
│ Records repeat, stub, copy, and topological features        │
└─────────────────────────────────────────────────────────────┘
```

---

## 2. Layer 1: Chromosome End State (`EndState`)

Evaluated once per chromosome arm end (e.g. `I:5p`, `IV:3p`). Captures the biological completeness of both somatic and germline lineage terminations.

| Enum Value | Description | Biological / Topological Role | Reference Benchmark Example |
| :--- | :--- | :--- | :--- |
| `end_complete` | Both somatic and germline arms capped | Authentic telomeric caps verified on both somatic terminus and germline arm. | Chr I 5p, Chr IV 3p (post-bridging) |
| `end_somatic_missing` | Germline cap present; somatic cap missing/stub | Authentic germline cap assembled; somatic diminution branch terminated in a low-support dropout stub or missing. | Chr V 3p |
| `end_germline_orphan` | Somatic cap present; germline cap absent or unbridged | Somatic telomeric cap verified; germline arm lacks an assembled or bridgeable terminal telomere. | Residual unbridged germline termini |
| `end_shared_exempt` | Single terminal cap shared by both lineages | End is exempt from PDE; somatic and germline lineages share the same distal telomere (e.g. distal to rDNA). | Chr II 5p |

---

## 3. Layer 2: Diminution Junction Classification (`JunctionClass`)

Classifies each Programmed DNA Elimination (PDE) junction or diminution break in the assembly graph.

| Enum Value | Description | Biological / Topological Role | Reference Benchmark Example |
| :--- | :--- | :--- | :--- |
| `pde_terminal` | Subtelomeric somatic diminution fork | High-coverage somatic stem unitig transitions into a somatic terminus and a lower-depth ($0.04\text{--}0.33\times C_{\text{core}}$) germline-retained sequence (GRS). | Chr I, II 3p, Chr III, IV, V, X |
| `pde_intercalary` | Internal chromosome fragmentation break | Internal chromosome diminution boundary where internal germline-retained material is excised during development. | Chr X Breaks 1 & 2 |

---

## 4. Layer 3: Interval Biological Quirks (`IntervalQuirk`)

Documents sequence-level and graph-level structural idiosyncrasies along chromosome intervals.

| Enum Value | Description | Biological / Topological Role | Reference Benchmark Example |
| :--- | :--- | :--- | :--- |
| `sat_gap_variable` | Subterminal variable satellite repeat array | Tandem repeat array spanning between core unitig and distal cap with copy-number disagreement across spanning ONT reads. | Chr IV 3p (~8.5 kb), Chr V 5p (~38 kb) |
| `sat_gap_fixed` | Fixed-length assembly/repeat gap | Repeat or assembly gap with deterministic spanning length across all reads. | Internal repeat transitions |
| `rdna_cassette` | Multi-copy ribosomal RNA tangle | Collapsed high-depth 18S-5.8S-28S repeating cluster with multiple entry/exit dovetails. | Chr II 5p |
| `uncollapsed_copy` | Low-depth alternate repeat copy/haplotype | Alternate repeat haplotype or paralog unitig surviving graph simplification at low read depth. | Chr II 5p satellite copy |
| `low_support_stub` | Low-support dropout stub | Terminal unitig failing read-support / depth criteria (e.g. 1–2 reads), representing somatic dropout. | Chr V 3p somatic stub |
| `bubble_rejoining` | Rejoined heterozygous/alternate bubble | Heterozygous bubble or alternate haplotype segment rejoining into the main stem unitig. | Heterozygous arm bubbles |
| `itr` | Internal telomeric repeat array | Canonical telomeric repeats trapped internally within a unitig or palindrome duplication, not serving as a cap. | PDE fork stem internal repeats |

---

## 5. Feature Identifier and Schema Conventions

### Canonical Feature Identifier Format

All records follow the canonical identifier syntax:
```
<feature_class>:<chrom>:<side>[:<index>]
```
Examples:
- `end_complete:I:5p`
- `pde_terminal:IV:3p`
- `pde_intercalary:X:mid:1`
- `sat_gap_variable:IV:3p:1`
- `rdna_cassette:II:5p:1`

In accordance with `BIO-006` and project anti-overfitting rules, **no unitig names or species identifiers may be hardcoded into classification logic**. Unitig IDs appear strictly in the evidence payload fields.

### Disambiguation of Low-Depth (1–4×) Unitigs

To distinguish authentic low-coverage distal germline caps (e.g. 1.0× with 7 spanning reads across 30 kb) from spurious dropout stubs (e.g. 1.0× with 2 reads), records incorporate:
- `read_count`: Number of distinct sequencing reads supporting the segment (from GFA `A` lines or read-association tags).
- `core_coverage_bp_median` ($C_{\text{core}}$): Median read depth across core somatic chromosome segments.
- `core_depth_ratio` ($D / C_{\text{core}}$): Normalized depth ratio relative to somatic core coverage.

---

## 6. Resolution Protocol for `sat_gap_variable`

When an assembly arm is interrupted by a subterminal variable-length satellite repeat array:

1. **Identification and In-Tool ONT Bridging**:
   - Detect spanning Oxford Nanopore (ONT) reads anchoring both the left unitig port and the right distal cap unitig port.
   - Extract the empirical length distribution across all spanning reads:
     - `spanning_read_count`
     - `length_min_bp`
     - `length_median_bp`
     - `length_max_bp`
     - `length_iqr_bp`
     - `repeat_monomer_len_bp`
     - `copy_number_range`

2. **Representative Sequence Selection (BIO-004)**:
   - Select a single representative spanning repeat sequence from the empirical read pool (typically matching the median read length) to bridge the interval.
   - **Zero `N` Gaps**: Do NOT collapse the repeat into a 100-bp `N` scaffolding gap. The finished assembly sequence must contain the authentic representative nucleotide sequence for downstream Tandem Repeats Finder (TRF) annotation.

3. **End-State Resolution**:
   - Once the distal germline cap is connected via the representative sequence, the chromosome end state resolves to **`end_complete`** (NOT `end_germline_orphan`).
   - `end_germline_orphan` is strictly reserved for germline caps that remain unbridged or unconnected even after bridging attempts.

4. **Quirk Ledger Record Emission**:
   - Emit a `sat_gap_variable` record in `quirks.tsv` documenting:
     - Feature ID (e.g. `sat_gap_variable:IV:3p:1`)
     - Bounding ports (`left_anchor_port`, `right_anchor_port`)
     - Read-length distribution metrics (`length_min_bp`, `length_median_bp`, `length_max_bp`, `length_iqr_bp`)
     - Repeat monomer details (`repeat_monomer_len_bp`, `copy_number_range`)
     - AGP metadata (`agp_component_type`, `agp_linkage_evidence`)
     - Representative sequence selection flag (`representative_sequence_selected = True`)
