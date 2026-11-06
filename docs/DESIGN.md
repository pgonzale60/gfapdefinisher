# gfapdefinisher: design brief (prototype)

Status: scaffold only. No modules or tests exist yet. This document is the handover of the architecture agreed in the design session. Home of the toolkit is undecided; it lives in its own repo for development.

## 1. Goal

A harness on top of hifiasm. It reads the hifiasm string graphs (`p_utg`, `r_utg`) and uses pre-computed biological annotations (USCOs, rDNA, telomeres, contaminant calls) to reach T2T germline chromosomes. It is not a de novo assembler. It must benefit from hifiasm's own graph refinement, so it never re-implements read walking.

It unifies two formerly separate efforts:
- chromosome-end resolution beyond somatic telomere junctions (PDE forks, rDNA, satellites);
- uncollapsed-allele / haplotype resolution in outcrossers (USCO + depth + topology).

Prototype scope: reach DF5120 (`nxOscSper1`) T2T automatically. DF5120 is a validation case, not a target. Nothing may be tuned to it.

## 2. Principles (non-negotiable)

1. **Every chromosome end is resolved.** An end that runs into rDNA or a satellite array is traversed to its distal telomere (legacy pdefinisher: rDNA precedence, satellite-linked telomere scaffolding). Never "accept" a naked rDNA end. rDNA end for DF5120 chr II is unsolved and will be solved during development.
2. **No arbitrary thresholds.** Fork classification must come from sequence composition (genuine terminal telomere vs internal telomeric repeats vs other), depth contrast and topology, not from lengths such as "stub < 10 kb". Any parameter that remains must be named, justified and read from one policy file.
3. **PDE breakpoints are inferred from the graph.** There is no external breakpoint BED. The fork is the evidence.
4. **No full read-set alignment is precomputed in production.** Round 1 inspects the graph, then extracts only relevant reads (hifiasm `A` lines, flank baiting) for targeted mapping. For prototype speed, fished alignments are precomputed (see section 5).
5. **Variable-size satellite gaps keep representation.** DF5120 chr IV and V 3' gaps are true variable-length tandem repeats (reads disagree on size). Do not emit a single N gap that erases the repeat. Emit enough repeat sequence for TRF to call a tandem repeat, and record the gap/variability in the AGP and audit.
6. **Lossless graph.** Preserve GFA overlap CIGARs, `rd:i:` and `A` lines. Do not convert to a blunt variation graph (odgi/vg).
7. **Phased subgenome paths are the ideal outcrosser output**, feasibility unknown. Design data structures so a walk can carry haplotype/subgenome labels; do not promise it yet.

## 3. Architecture

- **Graph core**: bidirected segments (two ends per unitig), `S`/`L`, CIGAR-aware overlap excision (exact `kM` and indel CIGARs).
- **Provenance store**: `A` lines streamed into an indexed on-disk table (SQLite/DuckDB), not Python heap. Query by segment, read id, haplotype tag.
- **Annotation layer**: pre-computed tables projected onto nodes (telomere tracts, rDNA, USCOs, contaminant class). Annotations are inputs, never computed inside the tool.
- **Classifiers**: fork/terminus classification; haplotype/bubble classification.
- **Walker + bridge engine**: path solving, targeted read fishing, abPOA consensus or audited repeat-preserving gap.
- **Emit**: FASTA + AGP v2.1 + audit (junction table, gap table, per-end evidence).

Optional compiled core (Rust/PyO3) was suggested by the "GFA String Graph Tool Evaluation" doc for scale. Not needed for DF5120 (~45 Mb); revisit only if larger genomes demand it.

## 4. Prior art in other repos: port vs avoid

Locations: `PDE_origin/scripts/assembly/hifiasm_t2t/`, `PDE_origin/scripts/analysis/usco_haplotype_toolkit.py`, `sanger_curation/nf-core-pdefinisher/`. Read the code before reusing; line ranges from the survey may have drifted.

**Reuse (adapt, add tests):**
- GFA `S`/`L` parsing, strand-aware edge storage (`graph.py`); overlap verification and merge (`overlap.py`).
- `A`-line / `rd:i` parsing and the depth-peak calibration for `asm.log` (hom/het peak quirk when `peak_het` is -1) in `usco_haplotype_toolkit.py`.
- Bidirected reachability with homologous-pair exclusion (bubble classes) from the USCO toolkit.
- Inverted telomere palindrome detection (`telomere.py`), as an artifact check.
- abPOA gap consensus and flank trimming (`gap_resolver.py`); AGP block accumulator and `strip_terminal_gaps` (`pde_generate_agp.py` in nf-core-pdefinisher).
- Terminus resolution hierarchy and the rDNA cassette concept from nf-core-pdefinisher, as the statement of intent for end resolution.

**Do not carry over:**
- `net_ext < 10000` somatic-stub rule and the hand-weighted candidate score (`extender.py`).
- GC window host filter (`contamination.py`); replaced by Tiara calls.
- Hardcoded DF5120 recipes and the embedded patch sequence (`recipes.py`).
- Fixed 2 kb terminal window and fixed 6-hop BFS limit; fixed minor-allele ratio 0.35.
- Any logic that treats an rDNA end as a legitimate terminus.

### 4a. nf-core-pdefinisher (the predecessor; read it first)

Repo: `/Users/pgonzale60/Documents/research/sanger_curation/nf-core-pdefinisher` (branch `dev`). It is the source of the project's intent: resolve absolutely every chromosome end. It works on scaffolded, semi-curated assemblies plus read BAMs, and models ends as open ports with valency conservation (`curationsynthesis/models.py`: `Port`). gfapdefinisher keeps that goal but works directly on the hifiasm graph instead of on scaffolds. How the two relate long term (replacement, wrapper, or a graph front end feeding it) is undecided until a working version exists.

Pieces that state the intended behaviour, to read and mirror at the logic level:
- `terminus_resolver.py`: the terminus hierarchy (rDNA array precedence, somatic boundary vs germline tip, unitig/read capping, satellite-linked distal telomere, uncapped subtelomere walk). Reimplement against graph evidence, not BAM clips.
- `rdnacassetteresolver/.../pde_rdna_cassette_resolver.py`: single-copy anchor + rDNA monomer + distal telomeric unitig, with a complexity gatekeeper (`DEFER_COMPLEX_RDNA_LOCUS`). Basis for the cassette model to be built here.
- `telomerewalker/.../teloscaffold/repeat_scaffolder.py` and `telocapping/unitig_capper.py`: satellite-linked telomere scaffolding and unitig capping with overlap excision.
- `targetedgapassembly/.../consensus_bridge.py`: spanning-read candidates and abPOA bridging.
- `curationagpgenerator/.../pde_generate_agp.py`: AGP v2.1 block accumulator.
- `assets/curation_policy.json`: the single source of truth for parameters ("policy parameters must not be hardcoded"). Copy this pattern. Caution: its values (e.g. `min_repeats` 15, `terminal_window_bp` 10000, `chimeric_terminus_max_dist_bp` 50000) are still fixed numbers and have not been shown to be non-arbitrary; each one adopted here needs a justification or a data-derived replacement.
- Its pipeline stages (reconnaissance, synthesis, resolution engines, surgery, evaluation) are the end-to-end flow gfapdefinisher must eventually be able to replace or feed.

## 5. Inputs for DF5120

Local: `data/benchmark/DF5120/` (gitignored; re-fetch from Farm22 if missing).
Farm22 hifiasm dir: `/lustre/scratch122/tol/teams/blaxter/users/pg17/nf_scratch/assemblies/Oscheius_sp_DF5120/results/nxOscSper1.pri/hifiasm/`
Farm22 precompute dir: `.../Oscheius_sp_DF5120/analysis/gfapdefinisher_precompute/`

| File | Source | Notes |
| --- | --- | --- |
| `asm.p_utg.gfa.gz`, `asm.r_utg.gfa.gz`, `asm.p_ctg.gfa.gz`, `asm.log` | hifiasm | core graph; `asm.p_utg.fa.gz` also on Farm22 |
| `tiara_unitigs.txt` | Tiara 1.0.3 on `asm.p_utg.fa` | deep-learning contaminant call; Tiara drops sequences < 3000 bp by default, so short unitigs may have no call. Check coverage before relying on it |
| `compleasm_unitigs_full_table.tsv` | compleasm 0.2.8, `chromadorea_odb12`, on `p_utg` | S 97.34%, D 1.51%. Only `p_utg`; `r_utg` USCOs not computed |
| `barrnap_rdna.gff`, `barrnap_rdna_hits.fa` | barrnap 0.9 `--kingdom euk` | 435 features; cassette model to be built during development |
| `tidk_telomeres_telomeric_repeat_windows.tsv` | tidk 0.2.7 `search TTAGGC`, 1 kb windows | 0.2.7 is the version in the cache; confirm whether a newer tidk was intended |
| `hifi_to_flanks.paf`, `ont_to_flanks.paf`, `test_flanks.fa`, `gap3_abpoa_input.fa`, `gap3_consensus.fa` | earlier hifiasm_t2t validation | fished against specific gap flanks of an earlier run; regenerate if flank definitions change |
| `nxOscSper1.1_GRS.fasta` | curated GRS loci | recovery audit set (13 loci) |
| `nxOscSper1.germline_T2T_validated.fa` | hifiasm_t2t checkpoint | **not ground truth**: 2 satellite gaps, uncapped rDNA end |

Raw reads (not staged): HiFi `.../tol_nematodes/Oscheius_sp_DF5120/genomic_data/nxOscSper1/pacbio/nxOscSper1_DF5120_hifi.hifi.fastq.gz`; ONT `.../genomic_data/nxOscSper1/ont/nxOscSper1.ont.combined.fastq.gz`.

Known gaps: no `r_utg` USCO table; tidk run was TTAGGC only; Tiara short-unitig coverage unchecked.

## 6. Validation (not tuning targets)

Reference metrics from the earlier toolkit, for comparison only: 97.15% of 828,758 bp GRS recovered, 13 loci, theoretical graph ceiling 97.64%. Success is reaching comparable or better results without hand-set recipes or thresholds, and with every end capped or traced to a distal telomere.

## 7. First milestone

1. `gfapdefinisher.graph`: lossless loader for `S`/`L`/`A` with tests (small synthetic GFA fixtures plus the DF5120 files).
2. `gfapdefinisher.annotations`: loaders for the Tiara, tidk, barrnap and compleasm tables, projected onto nodes.
3. A report listing, for every unitig end, its telomere/rDNA/fork evidence. No walking yet.

Then: fork classifier from composition + depth + topology, then walker, then bridging, then emit.
