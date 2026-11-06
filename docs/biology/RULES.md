# Biological Truth Rules (BIO Rules)

This document is owned strictly by the biological expert (`@pgonzale60`).
AI worker agents **MUST NOT** edit this document or invent biological rules.
All classifiers, decision engines, and termination logic in `gfapdefinisher` must cite one or more of these `BIO` IDs in their code, audit logs, and pull requests.

---

## BIO-001: Unconditional Terminal Completion
- **Statement**: Every chromosome arm must be traversed until it is terminated by an authentic telomeric cap. An end that enters ribosomal RNA (rDNA) arrays or satellite arrays is considered **unresolved** until the traversal continues through to the distal telomeric cap.
- **Rationale**: An rDNA or satellite cluster is an internal repeat, not a biological telomere. Accepting a naked rDNA end as an assembly terminus is biological nonsense and produces an incomplete genome.
- **Falsifiable Criterion**: Any assembly walk that terminates at a node with high rDNA hit density without an authentic telomeric repeat array is classified as `INCOMPLETE` / `UNRESOLVED`.

---

## BIO-002: Composition & Stoichiometry Over Arbitrary Length Cutoffs
- **Statement**: Fork and terminus decisions must be driven by physical sequence composition (canonical telomeric motifs, internal telomeric repeat sequences, ribosomal RNA subunits) and quantitative read depth contrast ($rd_{\text{grs}} \sim 0.04\text{--}0.40 \times rd_{\text{stem}}$), never by arbitrary sequence length cutoffs.
- **Rationale**: Arbitrary length cutoffs (e.g. `stub < 10 kb`) are fragile heuristics that fail across species with varying burst sizes and repeat lengths.
- **Falsifiable Criterion**: No code in `src/gfapdefinisher/` may compare node length against a static literal threshold to decide whether a branch is somatic or germline.

---

## BIO-003: Graph-Inferred Programmed DNA Elimination (PDE) Breakpoints
- **Statement**: Programmed DNA Elimination (PDE) boundaries are inferred directly from assembly graph topology and depth stoichiometry. There is no external breakpoint BED file provided in production.
- **Rationale**: The somatic vs. germline depth drop at the junction in the Hifiasm string graph is itself primary empirical evidence of elimination.
- **Falsifiable Criterion**: The assembly engine must identify somatic diminution boundaries without reading external genomic coordinate BED files.

---

## BIO-004: Repeat Representation Preservation for Variable Satellites
- **Statement**: Where variable-length tandem satellite arrays exhibit read-length disagreement across alleles (e.g., Chr IV and V 3' loci), the assembly must emit sufficient representative repeat sequence for Tandem Repeats Finder (TRF) detection rather than collapsing the array into a single 100-bp `N` gap.
- **Rationale**: Collapsing a complex multivalent repeat into a generic scaffolding gap obliterates biological satellite diversity and misleads downstream structural analyses.
- **Falsifiable Criterion**: Satellite-linked junctions must preserve representative repeat units in the compiled FASTA and record gap uncertainty in the AGP 2.1 component metadata.

---

## BIO-005: Lossless Overlap and Provenance Invariance
- **Statement**: String graph dovetail overlaps (including complex CIGAR strings with indels), coverage depths (`rd:i:`), sequence lengths (`LN:i:`), and embedded read placements (`A` lines) must be preserved without loss.
- **Rationale**: Down-converting OLC graphs to blunt variation graphs destroys unitig coordinates, read haplogroup provenance (`HG:A:p/m`), and excision accuracy.
- **Falsifiable Criterion**: Graph serialization and round-tripping must retain all non-standard Hifiasm tags and exact CIGAR operations.

---

## BIO-006: Held-Out Generalization Guard
- **Statement**: All heuristic rules and calibration metrics derived during development on `DF5120` must generalize and pass falsification against the three held-out test species: `JU1382`, `PS1017`, and `PS2068`.
- **Rationale**: Overfitting algorithms to DF5120 idiosyncrasies defeats the goal of a universal outcrosser and germline T2T finishing harness.
- **Falsifiable Criterion**: Any classifier rule that succeeds on `DF5120` but causes aberrant termination or misassembly on `JU1382`, `PS1017`, or `PS2068` is considered falsified and rejected.
