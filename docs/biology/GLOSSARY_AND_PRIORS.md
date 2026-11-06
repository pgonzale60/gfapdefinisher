# Biology Glossary and Anti-Prior Manifest

This document overrides any generic biological or computational assumptions embedded in LLM training weights.
When working on `gfapdefinisher`, **the definitions and rules below supersede all external knowledge**.

---

## 1. Project Glossary: What Words Mean in This Codebase

- **Telomeric Cap (Genuine Terminal Telomere)**:
  A run of canonical nematode telomeric repeats (e.g. `(TTAGGC)n` or reverse complement `(GCCTAA)n`) positioned at the physical end of a unitig, with outward-facing strand orientation and near-zero outgoing dovetail edges in the string graph.
- **Internal Telomeric Repeat (ITR)**:
  A telomeric repeat sequence that occurs internally within an arm or is trapped within an inverted palindrome duplication artifact. An ITR is **not** a chromosome cap.
- **Somatic Diminution Boundary (PDE Junction)**:
  A fork in the assembly graph where a high-coverage somatic stem unitig transitions into two alternative paths:
  1. A somatic truncated terminus (ending abruptly or with somatic telomere addition).
  2. A lower-coverage ($0.04\text{--}0.33\times$) germline-retained sequence (GRS) continuing the chromosome arm.
- **Germline-Retained Sequence (GRS)**:
  Chromosomal material present only in germ cells and deleted from somatic lineages during development. In Hifiasm graphs constructed from mixed tissues, GRS segments have reduced read depth proportional to the germline cell fraction ($0.04\text{--}0.33\times$ of core somatic coverage, $C_{\text{core}}$).
- **rDNA Cassette**:
  A repeating multi-copy array of 18S-5.8S-28S ribosomal RNA genes. In the graph, it forms a collapsed, high-depth tangle with multiple entry and exit dovetails.

---

## 2. Pretrained Priors That Are WRONG in This Project

AI models frequently default to common literature heuristics. The following priors are **explicitly forbidden** in `gfapdefinisher`:

| Common Pretrained Assumption | Reality in This Project (Strictly Enforced) |
| :--- | :--- |
| **"Short unitigs (<10 kb) are just noise or somatic stubs."** | **WRONG.** Micro-GRS segments, linker sequences, and genuine short unitigs can be <5 kb. Do not filter or categorize nodes by arbitrary length thresholds. |
| **"High GC% or low GC% indicates bacterial contamination."** | **WRONG.** GC content varies across nematode chromosome regions (e.g. AT-rich satellites). Contaminant classification is delegated to deep-learning models (Tiara), not global GC filters. |
| **"An assembly graph should be converted to a variation graph (vg/odgi)."** | **WRONG.** Blunt variation graphs destroy non-zero dovetail CIGAR overlaps and embedded read placements (`A` lines). We operate strictly on OLC GFA1 string graphs. |
| **"If an arm reaches an rDNA array, that is a legitimate terminus."** | **WRONG.** Chromosomes end in telomeres, not rDNA. An arm entering an rDNA cassette must be bridged across the cassette to its authentic distal telomere cap. |
| **"Variable satellite repeat differences can be collapsed into a 100-bp `N` gap."** | **WRONG.** Variable tandem repeats must retain representative repeat sequence for TRF calling, accompanied by AGP uncertainty annotations. |
| **"DF5120 unitig IDs can be hardcoded as special-case fallback rules."** | **WRONG.** DF5120 is a calibration benchmark, not a hardcoded target. No code may reference `utg000...` IDs or DF5120 names. |

---

## 3. Escalation Protocol for Biological Ambiguity

If an AI worker agent encounters a situation where:
1. A decision requires choosing a biological threshold not defined in `src/gfapdefinisher/policy/`.
2. A graph topology contradicts `BIO-001` through `BIO-006`.
3. An assumption is needed regarding chromosome arm architecture:

**The agent MUST STOP immediately, file or comment on a GitHub issue with the `bio-question` label, and await human expert input from `@pgonzale60`. The agent MUST NEVER improvise or guess.**
