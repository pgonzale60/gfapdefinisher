# PS2068 graph and telomeric-read refinement

The reviewed draft is **65,854,104 bp** across six chromosomes. Ten of the 12 physical tips pass the default canonical telomere scan. The subsequent open-port review accepts the short 14-repeat X 5′ tract on direct read evidence, bringing the end assignments to 11. III 5′ remains open. The draft corrects the rDNA orientation on II and uses a gapped representation of the array. It retains the frozen chromosome scaffold structure outside the documented terminal edits and rDNA block.

## Inputs and graph-first result

The frozen reference is `nxOscOnir1.2.primary.fa.gz` under `/lustre/scratch122/tol/teams/blaxter/users/pg17/freeze/2022-02-17/data/Oscheius/Oscheius_onirici/assembly/curated/`. It has six nuclear records totaling 65,783,231 bp. The supplied `asm.p_utg.gfa.gz`, `asm.r_utg.gfa.gz` and `asm.p_ctg.gfa.gz` come from the hifiasm 0.25.0-r726 run in `/lustre/scratch122/tol/teams/blaxter/users/pg17/nf_scratch/assemblies/Oscheius_onirici_PS2068/results/nxOscOnir1.pri/hifiasm/`. Its command used HiFi and `--ul` long-read guidance. The primary GFA SHA-256 equals the staged benchmark's SHA-256, `9ded9024ba4cef78e71e10528bbf87d43f4dc5f67410b6499059238219e5dab5`.

The current `gfapdefinisher finish` run used that primary GFA with the staged Tiara and Barrnap annotations. With FASTA and AGP emission and `--nominate-super`, it produced 4,125 FASTA records totaling 179,803,340 bp. The break-site summary lists two scaffolds with 114,923 bp of core and 47,084 bp of GRS. The II scaffold is rooted at the known II 3′ cut around II:10,439,111. The other starts from a unitig mapping at X:25,414–97,017, in frozen zero-based coordinates. The graph-first run yields local junction candidates and no chromosome-scale replacement for the frozen assembly.

Minimap2 2.29-r1283 aligned both unitig sets and the graph-derived contigs to the frozen chromosomes, and aligned the reference flanks back to the graphs. Individual full telomeric reads came from the supplied `end_clusters` run under the PS2068 `reads_only/nxOscOnir1` directory. Targeted reads came from `nxOscOnir1.2.hifi.sorted.bam`. The BAM has 94.79× mean mapped coverage across the reference. Both graph sets use the same input data. The draft uses one original molecule per terminal sequence and a base-exact 200 bp reference splice for the direct-read edits. Raw graph and read sequence are identical for the selected II 3′ extension. The X 5′ raw unitig also reproduces the selected HiFi overhang.

## Chromosome II rDNA and outward cap

The first 1,183,370 bp of frozen II contains 39 rDNA fragments and 39 200-N gaps before the core-side segment. The two complete or nearly complete monomers next to the core map in the reverse direction relative to the canonical monomer FASTA. Raw unitig `utg000013l` independently carries about 13.8 kb of rDNA in that direction and aligns from the rDNA entry into the next 168 kb of II. Within the gapped array, 20 fragments have the opposite orientation and 19 already agree with the core-side orientation. The draft reverse-complements the 20 discordant fragments in place. All 149 qualifying monomer matches in the corrected II audit have the same orientation. The original fragment order remains a scaffold representation.

`C000002` contains 108 rDNA-proximal telomeric molecules. The cap source `E00001173` (`m64016_201225_182119/100206578/ccs`) contributes an 886 bp telomeric tract and 4,704 bp of adjacent rDNA. The inward rDNA has the same orientation as the core-side and corrected array. This corresponds to an outward 5′ cap on II. A new 200-N join records uncertainty between that read fragment and the assembled array. The core-facing sequence is retained exactly.

## Terminal edits

The sizes below are changes relative to the frozen reference. The II 5′ addition includes its new 200-N join. The four terminal replacements preserve a unique reference anchor and use the observed read sequence through the cap. At I 5′ this closes the original 40-N terminal gap.

| End | Evidence family or source | Added bp | Retired reference bp | Result |
| --- | --- | ---: | ---: | --- |
| I 5′ | `C000034` | 8,778 | 483 | Outward cap from one of three exact splices. |
| I 3′ | `C000039` | 949 | 81 | Outward cap distal to the `C000010` elimination site. |
| II 5′ | `C000002` | 5,790 | 0 | Outward cap and reoriented array. |
| II 3′ | `C000014` and both unitig graphs | 7,353 | 0 | Exact cap sequence shared by the HiFi molecule and both graphs. |
| III 5′ | Full HiFi BAM | 1,393 | 0 | Single-read extension. `C000009` is at the internal elimination cut. |
| III 3′ | `C000025` | 7,943 | 0 | Outward cap on a complete selected molecule. |
| IV 5′ | `C000015` | 6,554 | 0 | Outward cap on a complete selected molecule. |
| IV 3′ | `C000032` | 2,365 | 140 | Outward cap beyond the `C000003` elimination site. |
| V 5′ | `C000017` | 9,082 | 837 | Outward cap beyond the `C000008` elimination site. |
| V 3′ | `C000006` | 10,842 | 0 | Outward cap. Both graphs also carry a longer telomeric version. |
| X 5′ | Full HiFi BAM and raw unitig `utg002074l` | 1,676 | 0 | A physical 14-repeat tract, below the 30-repeat default threshold. |
| X 3′ | `C000030` | 9,689 | 0 | Outward cap. The HiFi read ends at the tract and the longer graph continuation extends beyond it. |

The two `C000040` members illustrate why a family match alone leaves X 5′ short. `E00001388` reaches 7,809 bp into frozen X in the expected orientation and its extension matches raw `utg002074l` exactly. `E00001865` has a 12,073 bp telomeric tract but just 554 bp of inward sequence. That flank has a reverse alignment 4.5–5.0 kb inside frozen X and two placements within the raw unitig, so it cannot supply a uniquely anchored replacement for the short X 5′ cap.

## Reference conflicts and remaining gaps

The III unitig graphs place a 14.6 kb block from the frozen III 3′ arm in reverse orientation next to III 5′. Five reads span that candidate junction in the unitig panels. Four of the same reads also have a high-MAPQ split placement on IV 3′. Separately, 17 reads support the frozen III 3′ boundary, six with an exact 400 bp junction. The two read groups support competing paths in the terminal region. The reviewed draft retains the frozen III 3′ placement and uses its independent cap read. The 22 kb III 5′ graph proposal includes the disputed opposite-arm block. The 1,393 bp III 5′ extension comes directly from a HiFi read, and the targeted telomeric-read comparison supplies no anchored distal cap beyond it.

The draft retains 43 gaps of 200 Ns. Forty are in the reoriented II rDNA array. The other three are at the frozen I:6,135,625–6,135,825, II:2,952,558–2,952,758 and II:8,405,576–8,405,776, in zero-based half-open coordinates. None has an exact graph or single-HiFi-read bridge to both existing flanks in this audit.

## Verified outputs

The run directory is `/lustre/scratch122/tol/teams/blaxter/users/pg17/tol_nematodes/Oscheius_onirici/analysis/nxOscOnir1.2/ps2068_graph_refine.ICRF45/`. It contains the BGZF-compressed and indexed `PS2068.graph_read_refined.draft.fa.gz`, its per-end TSV, edit manifest, reproduction scripts and checksums. `II.rdna_corrected.json` records the 39 source intervals, strand changes and chosen cap. The baseline finisher logs and end calls are in `baseline_finisher/`.

## III and X 5′ open-port review

**X 5′ is accepted on direct evidence.** HiFi read `E00001388` anchors exactly across 7,809 bp of frozen X at MAPQ 60. Its 1,676 bp extension is identical to the raw graph, and its 14-unit telomeric tract ends one base from the physical draft tip. This is an end-specific curation call. The 30-repeat default scan threshold remains in effect.

**III 5′ remains open.** All cluster consensuses, representatives and full inward reads were compared against the 12 tips oriented from the outside toward the core. The verified 1,393 bp III 5′ addition has about 1.1 kb of related sequence at III 3′ (94.4% identity). It also has matches on IV and X, including a 1.39 kb reverse-oriented match near the X 5′ cap at 97.6% identity. The best `C000040` consensus hit to the III 5′ extension runs in reverse relative to the cap-to-core direction. No full archived member matches at least 200 bp of the verified III 5′ tip at ≥93% identity in the expected orientation.

The deeper 4,499 bp raw-graph segment before the competing III 3′ block supplies one forward-oriented hit in its outer 3.1 kb. This is the already assigned X 5′ read `E00001388`, beginning 3,430 bp inside its inward flank. It leaves the adjacent III 5′ cap placement unsupported. A separate HiFi overhang could supply 998 more bp at III 5′ by a local exact 200 bp splice, but that splice occurs three times in the draft and the read also split-aligns to IV. The assigned sequence ends with the previously verified 1,393 bp extension. A 200-N satellite join still requires a compatible outward-facing III 5′ cap fragment.

The updated `PS2068.open_port_review.ends.tsv` and `.json`, checksums and alignments are in the `open_ports.uQQR5k/` subdirectory of the Farm22 run. They refer to the same verified FASTA and preserve the distinction between 11 curated end assignments and ten caps meeting the default scan.

### III 5′ inversion hypothesis

If a distal cap is independently assigned to III 5′ and its inward-facing repeat runs opposite to the placed tip, an orientation switch within the satellite is a candidate explanation. A future review scaffold can represent that inferred junction with 200 Ns and keep the telomeric tract facing outward. The position and sequence of the switch are unresolved.

The unused long-telomere member `C000040`/`E00001865` carries a 12,073 bp tract and a 554 bp inward flank. The entire flank matches the sequence immediately after the assigned X 5′ tract exactly in the expected orientation. Its III 5′ alignment is reverse-oriented and partial, with 444 matching bases in 466 alignment columns. There is also a reverse internal X placement. This leaves another plausible source for the molecule and makes the III 5′ inversion conditional on identifying an independently assignable cap. `PS2068.III_5p.inversion_hypothesis.json` records this assessment alongside the open-port results. III 5′ remains open in the verified draft.

## Best placements for all telomeric families

The complete set of 45 telomeric families was ranked against the 11-cap review FASTA. Minimap2 aligned every full inward member and both family sequence summaries to the whole reviewed genome. Each molecule counts once per candidate end. End-specific HiFi anchors require MAPQ ≥20, ≥98% identity, ≥1,500 matching bases and a called telomeric tract within 25 bp of the sequenced end. A similarly strong competing hit to the same part of a molecule removes its unique-placement vote. The rankings prioritize end-specific molecules, then recorded cap-donor evidence, shorter high-quality anchors and compatible member counts, with consensus and representative matches as secondary evidence. The direction of the alignment identifies the candidate 5′ or 3′ port.

For tissue assignment, the inferred cap-to-flank boundary was compared within 500 bp to the placed cap and to the catalogue of major somatic cuts translated into draft coordinates. The predominant internal `C000001` site is used as an inferred I 5′ somatic boundary. More distal projections are tracked as germline candidates when the corresponding cut or draft cap provides a reference. Shared II 5′ rDNA assignments require rDNA-proximal reads, a match to the monomer and a match to the terminal array. Five rDNA families have reverse-only array hits and are flagged as weak candidates for that shared pool. Three families have only interior sequence matches, and three others have near-end sequence matches in the opposite orientation.

| End | Best somatic families | Best germline families | Shared rDNA |
| --- | --- | --- | --- |
| I 5′ | `C000001` (candidate) | `C000034` | |
| I 3′ | `C000010` | `C000039` | |
| II 5′ | | | `C000002` and ten smaller candidates |
| II 3′ | `C000004` | `C000014` | |
| III 5′ | `C000009` | | |
| III 3′ | `C000011` | `C000025` (accepted), `C000022` and `C000043` (candidates) | |
| IV 5′ | `C000013` | `C000015` | |
| IV 3′ | `C000003` | `C000032` | |
| V 5′ | `C000008` | `C000017` | |
| V 3′ | | `C000006`, `C000012` (candidate) | |
| X 5′ | `C000007` | `C000040` | |
| X 3′ | `C000005`, `C000024`, `C000042` (ambiguous placement) | `C000030` | |

III 5′ is the best end for the somatic family `C000009`, with 157 end-specific anchors. Every family with a sufficiently supported germline placement ranks elsewhere. `C000040` ranks at germline X 5′. `C000025` is accepted at germline III 3′ on four end-specific anchors. IV 3′ remains in the evidence as a competing sequence match. At V 3′, `C000012` has seven end-specific anchors and projects about 621 bp beyond the draft tip. Its candidate germline assignment is independent of any decision to choose a cap sequence.

II 5′ has **11 shared rDNA candidates**. Five singleton families match the array only in reverse orientation and retain the lower-confidence flag (`C000016`, `C000018`, `C000027`, `C000028`, and `C000033`). The reverse-only group may reflect a separate PDE-associated context at the rDNA end. The PDE interpretation remains a hypothesis until a locus-specific cleavage junction distinguishes an elimination cap from an alternative rDNA repeat placement.

The [45-family assignment table](PS2068_FAMILY_ASSIGNMENTS.tsv) includes the accepted placement and rankings with unresolved tissue or incompatible orientation. The full Farm22 output, including `family_end_rankings.tsv`, `family_assignments_by_end.tsv`, `family_curation_decisions.json`, the boundary translations and the alignment files, is in the `family_assignments.cmehfw/` subdirectory of the refinement run. The source FASTA matches its recorded checksum.
