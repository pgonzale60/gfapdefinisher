# gfapdefinisher design brief

The prototype has a working `finish` command with per-end evidence and FASTA and AGP outputs. DF5120 is the current validation case.

## 1. Goal

A harness on top of hifiasm. It reads the hifiasm string graphs (`p_utg`, `r_utg`) and uses pre-computed biological annotations (USCOs, rDNA, telomeres, contaminant calls) to reach T2T germline chromosomes. It is not a de novo assembler. It must benefit from hifiasm's own graph refinement, so it never re-implements read walking.

It unifies two formerly separate efforts:
- chromosome-end resolution beyond somatic telomere junctions (PDE forks, rDNA, satellites);
- uncollapsed-allele / haplotype resolution in outcrossers (USCO + depth + topology).

The current benchmark is DF5120 (`nxOscSper1`). Use it to evaluate germline end recovery without isolate-specific recipes.

## 2. Principles (non-negotiable)

1. **Germline end completion requires evidence of the distal telomere.** A telomere at a PDE cut can cap a somatic end while the germline tip remains missing. For DF5120 chromosome II, the core-side entry and telomere-side exit anchor the rDNA scaffold. The available reads are too short to resolve the array. Retain the assembler's rDNA copies in a consistent rDNA orientation and record uncertain joins as AGP gaps. The resulting scaffold represents the array between its flanks. Internal copy order, presence, copy count, and unknown gap lengths remain uncertain.
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
| `nxOscSper1.germline_T2T_validated.fa` | hifiasm_t2t checkpoint | Three satellite gaps and an uncapped rDNA end |

Raw reads (not staged): HiFi `.../tol_nematodes/Oscheius_sp_DF5120/genomic_data/nxOscSper1/pacbio/nxOscSper1_DF5120_hifi.hifi.fastq.gz`; ONT `.../genomic_data/nxOscSper1/ont/nxOscSper1.ont.combined.fastq.gz`.

Known gaps: no `r_utg` USCO table; tidk run was TTAGGC only; Tiara short-unitig coverage unchecked.

## 6. Validation (not tuning targets)

On the DF5120 baseline at `796d9f1`, `finish` used the `p_utg` graph, Tiara, Barrnap, and the staged ONT PAF and spanning reads. The six chromosome sequences carry telomeric tracts at both ends. The per-end evidence still has several open points.

- **IV 3′ and V 5′ in output coordinates.** The end ledger calls both core-unitig caps complete. Both ends have candidate germline continuations in the raw unitig graph (`r_utg`). Check the core-to-germline paths and the connections to distal telomeric sequence with ONT evidence. The output recovers only 0.2% and 15.2% of the corresponding curated terminal GRS loci, respectively.
- **II 5′.** Retain the oriented rDNA scaffold and its nine unknown-length AGP gaps. The run chooses the distal cap heuristically and marks the end `resolved=True`. Report support for the two flanks and the uncertainty of the array joins explicitly. The AGP places 397,594 bp before the core component, including the telomere-side component and gap placeholders. The cassette length fields report 8,311,179 bp because they include the core.
- **III 5′.** The satellite chain has two unknown-length gaps and zero spanning reads in the ledger. Its placement remains provisional.
- **Contamination (resolved).** The baseline emitted 157 Tiara-bacterial unitigs totaling 31,754,011 bp as standalone single-component scaffolds. None sat inside a chromosome. `finish` no longer seeds scaffolds from Tiara contaminants. It writes them to `excluded.fasta` with `excluded_by=tiara:<class>`, which leaves 15 scaffold records for DF5120. A resolution can still place a contaminant-called unitig inside a chromosome path. In PS1017, `utg000105l` (Tiara bacteria) is the germline branch of the `utg000016l` 5′ PDE fork: it has a 7,826 bp dovetail to the core and 15× depth against 50×. Its sister stub is also called bacterial. These calls suggest that Tiara misclassifies some eliminated sequence.

Using minimap2 `-x asm5` and merging query intervals at MAPQ ≥ 20, the six chromosome sequences cover 731,792 of 828,758 curated GRS bases (88.3%) across 13 loci. The same calculation gives 545,904 bases (65.9%) on the checkpoint assembly. The earlier toolkit's 97.15% used an undocumented calculation.

## 7. Benchmark order and Hi-C evaluation

Finish the DF5120 investigation first. Use high-coverage JU75 as the next case. Its close relationship to CEW1 and the expectation of end-only DNA elimination give a hypothesis to test when auditing chromosome ends and checking for internal GRS.

The [DF5120 graph refinement audit](DF5120_GRAPH_REFINEMENT.md) compares the existing chromosome FASTAs against both hifiasm unitig graphs. The [DF5120 terminal evidence audit](DF5120_TERMINAL_EVIDENCE.md) traces the historical curation across all 12 ends and explains which placements the staged checkpoint and regression test miss.

Keep graph-guided refinement of existing Hi-C-linked chromosome assemblies open as a design option. Start the comparison from the reviewed chromosome scaffolds and use the graph and read evidence to identify specific candidate corrections. In PS2068, account for retained allelic variants and preserve prior allele choices unless the evidence supports changing them. Compare the corrections and their support with the graph-first workflow before choosing a direction.
