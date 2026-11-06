# gfapdefinisher

**Biologically Informed Hifiasm String Graph Harness for Automated Germline T2T & Outcrosser Assembly Finishing**

`gfapdefinisher` operates directly on Hifiasm GFA1 assembly graphs (`asm.p_utg.gfa.gz`, `asm.r_utg.gfa.gz`), preserving overlap CIGAR strings, coverage depths (`rd:i:`), and embedded read placements (`A` lines) without destroying string graph coordinates.

## Core Design Principles

1. **Unconditional T2T Completion**: Every chromosome arm is traversed past programmed DNA elimination (PDE) junctions, through repeat arrays, and across rDNA cassettes until terminated with authentic telomeric caps.
2. **Composition & Stoichiometry over Arbitrary Cutoffs**: Fork decisions are driven by physical sequence content (telomeric repeats, internal telomeric repeat sequences, microhomology) and quantitative somatic vs germline depth contrast ($rd_{\text{grs}} \sim 0.20\text{--}0.40 \times rd_{\text{stem}}$), completely discarding arbitrary length cutoffs.
3. **Decoupled Architecture**:
   - **In-Memory OLC String Graph**: Bidirected dovetail edges, coverage depth, and biological interval overlays.
   - **Out-of-Core Read Index**: Streaming parser for millions of `A` lines into an indexed SQLite/DuckDB cache, queried on-demand for phase consistency and targeted gap fishing without Python heap bloat.
   - **Pre-computed Biological Layers**: Deep learning contaminant masks, single-copy orthologs (compleasm / BUSCO), and ribosomal RNA arrays.
   - **Exact CIGAR Overlap Excision & AGP 2.1 Surgery**: Emits verified FASTA alongside authoritative INSDC AGP v2.1 component maps.
