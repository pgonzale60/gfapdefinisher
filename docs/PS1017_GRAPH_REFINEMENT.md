# PS1017 graph and telomeric-read refinement

PS1017 is *Oscheius dolichura*. The reference-guided run starts from the six curated nuclear chromosomes in `nxOscDolc1.1.primary.fa.gz`. The evidence for retained junctions includes primary and raw hifiasm unitigs, primary contigs, original HiFi molecules, and a targeted check of the archived ONT BAM. The verified draft has **69,976,127 bp**, **seven assigned distal caps**, and **one unresolved 200-N join on X**.

## Inputs and graph-first baseline

The frozen reference is at `/lustre/scratch122/tol/teams/blaxter/users/pg17/freeze/2022-02-17/data/Oscheius/Oscheius_dolichura/assembly/curated/nxOscDolc1.1.primary.fa.gz`. The six nuclear records total 69,661,230 bp. The FASTA contains the six nuclear chromosomes. The eight frozen unlocalized records served as additional alignment targets in the family audit.

The hifiasm 0.25.0-r726 outputs are in `/lustre/scratch122/tol/teams/blaxter/users/pg17/nf_scratch/assemblies/Oscheius_dolichura_PS1017/results/nxOscDolc1.pri/hifiasm/`. The hifiasm command used HiFi and `--ul` guidance. The primary unitig GFA matches the staged benchmark by SHA-256, `79155d323d903b64761fe0c18633c090e581e5c31551b2d1558efa7d9d2154bb`.

With the staged Tiara and Barrnap annotations, `gfapdefinisher finish --nominate-super` produced 236 FASTA records totaling 84,063,386 bp. The break-site summary has 16 SUPER scaffolds and one smaller scaffold, with 26,298,931 bp classified as core and 7,360,677 bp as GRS. The run reported 13 caps and 21 break sites among those scaffolds. Several proposed GRS blocks had depth ratios of 0.82–1.35. The previously measured average coverage is 63.7× for core sequence and 10× for GRS. These results led to refinement of the frozen six-chromosome assembly.

The alignments use minimap2 2.29-r1283. The run used HiFi BAM `analysis/nxOscDolc1.1/read_mapping/nxOscDolc1.1.hifi.sorted.bam` and the archive of 56 telomeric families in `analysis/reads_only/nxOscDolc1/telomeric_read_fishing/end_clusters/`, under the Farm22 `tol_nematodes/Oscheius_dolichura` directory.

## Gap and graph conflicts

The frozen assembly has three internal 200-N joins. At I:1,930,383–1,930,583, the replacement is a 51 bp allele found in all three graph exports and 60 distinct HiFi molecules, with two exact 125 bp flanks. Other reads carry longer satellite alleles. The raw graph also contains an inversion branch near this gap, with one spanning HiFi molecule in the tested panel. The selected allele follows the frozen chromosome orientation.

The primary contig `ptg000004l` joins an I block to V in reverse orientation. The contested frozen boundaries are near I:10,627,634 and V:3,069,710. Targeted checks found 32 reads continuously covering the I boundary and 18 covering the V boundary, with zero reads spanning the proposed cross-chromosome junction in that panel. The draft keeps the two frozen chromosome structures.

See `gap_exact_spans.tsv`, `gap_local_alternatives.tsv`, `split_placement_candidates.tsv`, `panel_read_support.tsv` and `frozen_structural_support.tsv` for the underlying counts. Coordinates in this report refer to the frozen reference and are zero-based unless specified.

## Terminal curation

| End | Result | Change in bp | Evidence |
| --- | --- | ---: | --- |
| I 5′ | Germline cap | +14,391 | The contig and oriented unitig paths support the continuation. An original `C000025` read supplies the distal tract. |
| I 3′ | Germline cap | +38,805 | Primary unitig `utg000049l` and `C000043`. Six HiFi reads cross the frozen tip and 15 cross into sequence matching `III_unloc_1`. A longer contig route has distinct `C000028` support and remains an alternate. |
| II 5′ | Open | 0 | Both oriented unitig routes stop at the frozen tip. The major `C000005` family starts about 67 kb inward, at the known elimination cut. |
| II 3′ | Germline cap | +4,017 | One `C000037` molecule spans the original tip and the cap, with an exact anchor shared with the graph. |
| III 5′ | Extended, open | +14,134 | Identical sequence across raw and primary unitigs and the primary contig. Read `m64016_201225_182119/72747416/ccs` supports most of the extension (11,874 bp beyond the old tip). |
| III 3′ | Germline cap | +9,868 | One `C000047` molecule spans the original tip and the cap. |
| IV 5′ | Germline cap | +1,449 | A `C000042` molecule reaches the terminal tract from a unique exact anchor. The edit replaces the outer 2,979 bp with a 4,428 bp observed allele. |
| IV 3′ | Outward rDNA cap | +2,973 | A distal `C000019` molecule supplies a continuous 7,944 bp terminal segment with a unique exact splice replacing the outer 4,971 bp. |
| V 5′ | Germline cap | +213,191 | Primary unitig `utg000042l` places part of `I_unloc_1` on V. Two original HiFi reads cross the frozen tip and three cross into the unlocalized block. A `C000039` read named in the GFA path supplies the outward tract. |
| V 3′ | Extended, open | +6,877 | One HiFi molecule spans the original tip into the repeat. A longer contig route traverses repeat copies and `V_unloc_2`. Its full repeat span and cap placement are uncertain. |
| X 5′ | Extended, open | +796 | Both unitig graphs replace 9 bp with 805 bp. The terminal repeat has competing placements and `C000001` has a stronger reverse match. |
| X 3′ | Extended, open | +4,235 | One HiFi molecule exactly spans the continuation shared by both unitig graphs. The longer primary-contig path runs into repeated sequence with conflicting placements. |

### IV 3′ rDNA

The canonical 7,584 bp monomer aligns to the core-facing array in forward orientation. The outward splice uses `E00001100` from `C000019`. Thirty-eight independently named `C000019` molecules agree with the known inward breakpoint at IV:12,459,187 in the frozen assembly, and eight project to the distal rDNA beyond the frozen tip. This supports both positions within the same family. Repeated copies limit the assignment of individual reads to a particular rDNA copy. The selected outward sequence is continuous and comes from an observed molecule. No rDNA gap was added.

The four rDNA-associated families are `C000009`, `C000019`, `C000052`, and `C000053`. The family table records the placement evidence, including competing loci and the direction of monomer matches.

## Archived ONT follow-up

The indexed ONT BAM is `/lustre/scratch122/tol/teams/blaxter/users/pg17/freeze/2022-02-17/data/Oscheius/Oscheius_dolichura/analysis/coverage/minimap/ont/nxOscDolc1.1.bam`. The targeted survey recovered 1,430 complete read sequences by name from indexed gap, terminal and unlocalized intervals.

- **III gap closed.** Three ONT reads in the targeted alignments reach both flanks, with MAPQ 60 on each side. The graph flank alignments terminate on separate raw unitigs and the ONT sequence aligns to multiple copies of the repeat. Read `04f16059-21c5-4fc8-a75e-5ea726a5d6c0` supplies 4,520 bp of continuous unpolished sequence at exact 42 bp and 39 bp splice anchors. The edit replaces frozen III:4,641,925–4,642,135, including ten bases just before the original gap. This gives a net addition of 4,310 bp. The other spanning molecules support alternative lengths. See `ONT.gap_exact_splices.tsv` and `ONT.gap_alleles.{p_utg,r_utg}.paf`.
- **X gap remains.** The original X:3,833,531–3,833,731 retains its 200 Ns. The 4 kb flanks terminate on separate raw unitigs, and the targeted HiFi and ONT alignments contain no molecule that bridges both. The gap is at X:3,834,327–3,834,527 in the review FASTA.
- **X 5′ remains open.** Two ONT molecules carry a terminal tract and extensive sequence related to the X 5′ repeat. Whole-assembly and raw-graph alignments place several blocks from each molecule at overlapping X positions in both orientations. The tract-adjacent block of `e10da57f-450f-4578-81e8-6093d1e4fa96` also aligns for 7.4 kb near the accepted I 5′ cap at MAPQ 36. The other read, `a95b8e0e-37c5-43cf-8926-5df12b73245d`, is a possible X 5′ cap candidate, with competing placements through the repeat. This evidence leaves the repeat junction ambiguous. See `ONT.X5.whole.paf`, `ONT.X5.raw_graph.paf` and `ONT.X5.overlaps.paf`.
- **The remaining open ends are II 5′, III 5′, V 3′ and X 3′.** The targeted ONT data did not establish an outward cap for these ends or span the full V 3′ repeat path into `V_unloc_2`.

## Family assignments and outputs

All 56 archived families were ranked against the final reviewed assembly and the frozen unlocalized records. Fifteen best assignments match known internal somatic cut sites. Twelve match end-associated somatic sites, and `C000019` has evidence for both positions at the IV 3′ rDNA array. The table gives best placements and confidence for families with limited discrimination. [PS1017 family assignments](PS1017_FAMILY_ASSIGNMENTS.tsv)

Farm22 run directory:

```text
/lustre/scratch122/tol/teams/blaxter/users/pg17/tol_nematodes/Oscheius_dolichura/analysis/nxOscDolc1.1/ps1017_graph_refine.MOnhff/
```

The final deliverables are `PS1017.review.fa.gz` with `.fai` and `.gzi` indexes, `PS1017.review.ends.tsv`, the coordinate and sequence-source manifest `PS1017.review.json`, and `PS1017.family_assignments.tsv` with full alternatives in `PS1017.family_rankings.tsv`. The checksum manifest covers the FASTA, both indexes and the reports. The previous HiFi-only draft is archived under `hifi_only_review/`. The run holds the audit scripts, exact splice evidence and build policy used to reproduce the final draft.

## Open-port audit of the II and III 5′ short satellites

The audit used the outer 23 bp repeat on II and the outer 20 bp repeat on III. The marker screen compared oriented chromosome and unlocalized ends, all three graph sequence sets, the archived telomeric reads, and primary sequences from the original HiFi and ONT BAMs. Candidate molecules were aligned back to the reviewed assembly and raw graph. Results are in `ii_iii_open_ports/` under the Farm22 run directory.

**II 5′ connects to `II_unloc_3`.** Its 3′ edge shares 89 of 124 distinct 21-mers from the II tip in the orientation expected for placing `II_unloc_3` before II. Four independent ONT reads anchor on both sides of the join at MAPQ 60 and span 10,728–10,991 bp between the two flank alignments. All four yield exact splice anchors. Read `6a278710-be23-46f4-a731-97ce7010c03d` supplies a 10,983 bp continuous interval with exact 104 bp and 39 bp anchors. The proposed sequence joins `II_unloc_3[0:39200]`, the ONT read interval `[18211:29194]`, and reviewed `II[4198:]`. It adds 45,985 bp to II. The bridge is directly sequenced and requires no Ns. The primary and raw unitigs account for portions of the two flanks, with no single graph sequence spanning the same interval. The ONT interval is unpolished.

The 90 kb bridge read `1338a409-35c5-4261-9b6e-7aa5bde734bd` extends beyond `II_unloc_3` and aligns in reverse orientation near the first 20.7 kb of `II_unloc_1`. This supports a farther route into that fragment, with its 3′ edge facing outward. The route lacks a verified outward telomere. `C000034` has competing placements in the shared sequence of `II_unloc_1` and `II_unloc_4` and offers no unique outward anchor at the predicted distal edge. The exact II splice is recorded as a proposal in `II23.exact_splices.tsv` and `II23.ONT_patch_alleles.fa`.

**III 5′ still has competing paths.** Primary unitig `utg000068l` and contig `ptg000020l` describe the same roughly 203 kb capped fragment. It contains `V_unloc_1`, and family `C000038` has reads that align across the adjacent cap flank. Its inward-facing edge carries the 20 bp repeat in a compatible orientation. The proposed V 3′ graph path reaches a homologous open repeat through `V_unloc_2`. Raw unitig `utg000162l` matches the last 33 kb of that V path at over 99.9% identity. The capped fragment shares 352 distinct oriented 21-mers with the repeat-facing V path and 269 with the current III 5′ port, using the comparison windows in `oriented_marker_comparison.tsv`. These findings make V 3′ a plausible competing placement for the cap. The targeted molecules supply no unique bridge from the capped fragment to III 5′. The strongest apparently telomeric HiFi reads directly aligned to the III 5′ arm have their telomeric-like signal on the inward side. Other screen hits localize to already-assigned ends or unlocalized repeated sequence.

## II 5′ extension through both unlocalized fragments

The revised FASTA in `ii5_chain_extension/PS1017.II5_extended.review.fa.gz` extends chromosome II by **128,537 bp** relative to `PS1017.review.fa.gz`. The indexed six-chromosome assembly totals **70,104,664 bp**. II 5′ is still open. This revision adds no gap to II.

The extension uses `II_unloc_3` in forward orientation and `II_unloc_1` in reverse orientation. The 90 kb ONT read `1338a409-35c5-4261-9b6e-7aa5bde734bd` joins their two satellite flanks. Its interval `[18070:20647]` supplies 2,577 bp between exact 59 bp and 54 bp anchors. The sequence starts with the first 76,415 bp of reverse-complemented `II_unloc_1`, that read interval, and `II_unloc_3[4373:39200]`. The previously validated 10,983 bp ONT splice links this sequence across the 23 bp satellite to II. These two splices add 120,604 bp. The replacement intervals are unpolished ONT sequences from individual observed molecules.

Beyond the exposed 3′ edge of the original `II_unloc_1`, primary unitig `utg000244l+` and raw unitig `utg000402l-` carry an identical **7,933 bp** continuation. Four HiFi reads cross the former tip, and two ONT reads cover the graph extension. This extension is incorporated in the new FASTA. Its outward edge remains in mixed satellite DNA with local 191 bp periodicity and longer repeat patterns around 725 and 916 bp. The graph sequence has no terminal tract.

ONT read `d35abc4d-e552-496a-a751-7bd7c730b0cc` aligns across the graph continuation and reaches about 8.5 kb farther. The candidate in `II5_graph_and_ONT_candidate.fa` retains 7,858 bp of graph sequence and adds 8,534 bp from that read at a 36 bp exact splice. The distal part matches several repeat placements in `II_unloc_2`. The primary graph also proposes a farther route containing `II_unloc_4` and part of `II_unloc_2` in `utg000063l`. The alternate path and ONT overhang remain candidates because their placement through the repeats is ambiguous. The evidence and sequence sources are in `II5_chain.draft.json`, `sequence_support.tsv`, `outward.far_route.paf` and `outward.whole.paf`.
