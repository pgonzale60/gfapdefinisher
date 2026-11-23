# JU75 reference and graph refinement

## Reference baseline

The registered Farm22 `nxOscSpeh2.1.primary.bgz.fa.gz` contains six nuclear chromosomes, one mitochondrial record and 56 unplaced records. The staged reference contains 59,562,621 bases. The existing HiFi coverage summary gives mean nuclear chromosome depths of 410–453×. This is depth against the frozen reference. A separate documented hifiasm attempt used a 100× capped input.

| Chromosome | Reference length | Internal N-runs | Catalogued PDE junctions |
| --- | ---: | ---: | --- |
| I | 8,699,894 | 1 | Both ends |
| II | 8,404,209 | 4 | 5′ |
| III | 9,750,874 | 0 | Both ends |
| IV | 10,788,826 | 4 | Both ends |
| V | 8,944,530 | 1 | 5′ |
| X | 10,882,673 | 2 | 3′ |

Eleven gaps are 200 Ns and one on IV is 23 Ns. The nine sites in the frozen PDE catalog are all terminal. The other chromosome ends need independent evidence before assigning a breakpoint count. The native telomere scan found no tracts meeting its terminal-cap criteria on this FASTA. II 5′ has a substantial telomeric tract starting 1,134 bp from the physical tip.

## Annotation for anchor-based curation

The frozen reference already has eukaryotic 18S and 28S annotations near the start of X. A 6,568 bp JU75 rDNA monomer aligns over its full length to three tandem copies spanning X:15,429–35,136 in 1-based coordinates. The monomer also hits the unplaced X-associated records. The later curation dossier places the array at X 3′ in its described assembly. Use the frozen orientation for any edits to `nxOscSpeh2.1` and identify the rDNA core-side entrance and the telomeric flank independently.

Start with the GFA links, overlap sequence, depth and terminal telomere scans at reference-anchored autosomal gaps and ends. Two unique chromosome flanks on a graph path provide the strongest placement check for an internal patch. For terminal paths, compare the distal branch to any somatic cap, and check overlaps and competing paths before extending into low-depth sequence. Constrain output to placed chromosome components. This contains the risk from disconnected contaminants. Ambiguous or weakly anchored extensions still need a targeted sequence or taxonomic check. Track 45S rDNA on X because placement there depends on distinguishing rDNA copies, intervening sequence and the distal telomeric component.

## Graph availability

The documented Farm22 hifiasm output directory for JU75 contains `asm.log` only. Five recorded attempts exited at a memory or runtime limit, including the 100×-capped Hi-C-integrated attempt (job 103254, memory limit). The original Tree of Life archive has the curated FASTA, a decontaminated hifiasm draft FASTA and curation metadata. Neither archive exposed a GFA. The PacBio CCS read BAMs are archived. A recovered or newly generated HiFi graph is needed to verify repairs or terminal paths.

## Read composition before a new assembly

In job 103254, hifiasm estimated a homozygous peak of 98× and indexed 766,726,844 minimizer positions. LSF killed the run after 4 hours 49 minutes at 64,004 MB against a 64,000 MB limit. The original assembled draft's FCS-GX action report has zero actionable contamination calls. Some taxonomy rows carry small bacterial hits, including a 455 bp prokaryotic match on a 17,316 bp repeat record. The historic screen establishes little about the unassembled or unmapped portion of the later read input.

Two aligned BAMs exist for the same two trimmed HiFi libraries. The copy under `read_mapping/` has zero unmapped reads because the mapping pipeline used `samtools view -F 4`. The 13 GB BAM under `analysis/nxOscSpeh2.1/coverage/nxOscSpeh2_JU75_hifi.sorted.bam` retained **2,032,421 unmapped reads** and 2,523,086 primary mapped reads. Its recorded somatic-core mean is **444.96×**. The two PacBio movie names in its minimap2 header are `m64230e_220603_030004` and `m64230e_230125_115038`.

The old 100× filter sampled chromosome-core reads and retained **2,032,112 unmapped records** after separate rDNA and mitochondrial filtering, about **78%** of the 2,610,681 selected records. It excluded secondary alignments from selection, while supplementary records could still affect classification. The capped FASTA index records 2,604,321 sequences totaling **36.57 Gb**, equivalent to roughly **636×** the 57.47 Mb nuclear reference length in input bases. This is an input-volume ratio. The mapped somatic-core depth was capped separately.

Reuse the full coverage BAM for read partitioning. Screen its unmapped reads first for positive bacterial evidence, using read-scale competitive coverage or taxonomic assignments. Use all host reference records and the JU75 45S rDNA monomer when recruiting host reads. Track read counts and bases for host, confirmed foreign and unclassified bins, and retain the unclassified set separately. Preserve chromosome-end and breakpoint-related reads and handle high-copy rDNA and mitochondrial reads separately. Pick the host-core sampling depth after quantifying those bins, then build the HiFi-only primary and raw graphs in a fresh directory. [Hifiasm coverage and memory options](https://hifiasm.readthedocs.io/en/latest/parameter-reference.html), [FCS-GX methodology](https://pmc.ncbi.nlm.nih.gov/articles/PMC10898089/)

## HiFi graph refinement draft

The HiFi-only hifiasm build used 3.26 Gb from the coverage BAM. Tiara 1.0.3 flagged 848,417 of 2,032,421 unmapped reads with two bacterial scores of at least 0.95. The unmapped rescue set was retained separately. Hifiasm 0.25.0-r726 produced a 1,009-unitig primary graph and a 1,021-unitig raw graph.

Minimap2 2.29-r1283 aligned both unitig graphs and the graph-derived contigs to `nxOscSpeh2.1`, and aligned the chromosome termini and gap flanks back to each graph. The draft uses graph sequence with 160 bp of base-exact reference context at each splice. The HiFi checks used targeted reads from the full coverage BAM. Terminal span counts below require MAPQ ≥20, alignment identity ≥98%, 500 bp of reference-side anchor, and an alignment reaching within 100 bp of the graph tip. The two unitig graphs share the same underlying HiFi data.

The six-chromosome draft is **57,515,843 bp**. Ten of the 12 internal gaps closed at graph-supported overlaps of 1–17 bp. Nine have exact sequence support in both unitig graphs. At II gap 2, a 10 bp discrepancy on the right reference flank was removed before merging a 9 bp overlap. There were 143 HiFi reads that matched the revised junction exactly across 100 bp. The draft adds 46,736 bp of graph sequence at seven chromosome ends.

| Reference end | Change | HiFi spans | Result |
| --- | ---: | ---: | --- |
| I 5′ | +1,268 bp | 16 | Short raw-unitig extension into the subtelomere. The end remains open. |
| I 3′ | 0 | | The short graph cap lies near the catalogued PDE cut. Distal germline placement remains open. |
| II 5′ | 0 | | The reference has a telomeric tract 1,134 bp from the tip. The graph alternative would replace distal reference sequence. |
| II 3′ | +7,882 bp | 6 | Telomeric cap added. |
| III 5′ | +7,062 bp | 1 | Telomeric cap added. The longer raw alternative has no full span in the targeted alignments. |
| III 3′ | +1,935 bp | 2 | Telomeric cap added from an exact two-unitig raw path. |
| IV 5′ | +5,139 bp | 4 | Telomeric cap added distal to a separate unresolved gap. |
| IV 3′ | 0 | | The graph cap would replace distal reference sequence. A competing alignment places sequence from the reference tip on III. |
| V 5′ | +5,149 bp | 1 | Telomeric cap added. The spanning read maps to V at reference position 3,257 with MAPQ 60 in the source BAM. |
| V 3′ | +18,301 bp | 20 | Telomeric cap added from the shorter raw path. The 41,450 bp contig route has no full HiFi span and includes shared repeats. |
| X 5′ | 0 | | The graphs supply an rDNA-associated cap candidate. No targeted read from the placed reference anchor reaches the telomeric tract for either candidate. |
| X 3′ | 0 | | The forward graph cap starts inside the reference tip. The end-reaching match reverses orientation. |

Two gaps remain, both 200 Ns. The I gap is at `6,931,815–6,932,015` in the draft and the IV gap is at `26,546–26,746`, in zero-based half-open coordinates. Neither graph has an exact bounded path between the I gap boundaries. Both graphs spell a 26,051 bp path at IV 5′. The targeted HiFi reads support each side separately, with no read spanning that interval and both flanks. The new IV cap therefore establishes the distal segment and does not establish a continuous path from the core. The six new telomeric tracts record sequence evidence at the physical tips. Germline status and continuity from the chromosome core still require per-end review.

The Farm22 run is under `/lustre/scratch122/tol/teams/blaxter/users/pg17/tol_nematodes/Oscheius_sp_JU75/analysis/nxOscSpeh2.1/ju75_hifi_rescreen.dsGoEq/refine.7W80zD/`. `JU75.graph_refined.draft.fa.gz` is BGZF-compressed and indexed. `JU75.graph_refined.ledger.tsv` records all 24 features. `JU75.graph_refined.edits.json` records source paths, oriented graph coordinates and insert checksums. The run also contains the reviewed choices, the patch application script and the read-evidence tables.

## Telomeric read clusters

The cluster members from the Farm22 `telomeric_read_fishing/end_clusters/` run were aligned as individual full inward reads to the six-chromosome draft, frozen unplaced sequences, and the primary and raw graphs. The uniquely anchored placements use a telomere within 15 bp of the sequenced tip, MAPQ ≥20, at least 1,500 aligned matching bases, and alignment identity ≥97%. The 335 sequence families are clustering results. Sixty-one have at least one strong nuclear placement within 250 kb of a draft end.

Three alternate caps have both a unique, exact 500 bp draft splice and original telomeric reads covering that splice on the raw graph. Each proposed splice replaces existing distal draft sequence. The I 3′ cluster `C000173`, with 371 reads matching the graph cap, anchors about 2.4 kb inside the draft tip near the catalogued elimination cut. Neither I end has a supported distal connection in this read set.

The dominant `C000173` junction leaves 2,411 bp distal to it in the frozen assembly, at I:8,697,484–8,699,894 (1-based, inclusive). Earl Grey `nxOscSpeh2.1.primary.filteredRepeats.gff3` annotates an 820 bp `Unknown` repeat (`RND-1_FAMILY-62`) at 8,698,255–8,699,074 and a 762 bp `LTR/DIRS` fragment (`RND-1_FAMILY-112`) at 8,699,122–8,699,883. The satellite tracks have no calls in this interval. The physical tip carries just 18 bp in a `TTAGGC` phase, equivalent to three repeat units.

| End | Families | Spanning HiFi reads | Draft bases replaced | Graph bases added |
| --- | --- | ---: | ---: | ---: |
| II 5′ | `C000169` | 126 | 2,668 | 7,270 |
| IV 3′ | `C000181`, `C000319` | 12 | 11,622 | 9,749 |
| X 3′ | `C000180`, `C000262` | 54 | 20,873 | 7,500 |

At X 5′, two full `C000161` members map through the distal raw unitig `utg001015l+` and into `X_unloc_25`. One reaches all the way through the selected set of four complete distal monomers. Two raw walks from placed X lead through rDNA to caps matching the same family. The walks share a 21,979 bp nominal overlap that differs at 10,103 positions, and the read-backed `utg001015l` has no graph links. These observations support a provisional X 5′ placement with the exact join and array copy count open. The review candidate replaces the first 2,298 bp of the draft with 35,376 bp from `utg001015l+` and a 200-N gap at a common monomer boundary. The distal fragment and the five complete monomers retained on placed X run in the same orientation. Their combined nine complete copies are an assembly representation.

Neither internal gap gains a chromosome-anchored spanning read. Two `C000185` reads align across nearly the entire rejected 26,051 bp IV insert, but both map at MAPQ 60 across more than 30 kb of V. That placement rules out using those reads to close the IV gap. A third read matching part of the insert also maps uniquely to V.

The separate `JU75.cluster_connected.review.fa.gz` applies the four candidate ends and has ten physical telomeric caps. Both I ends remain open. Three 200-N gaps remain at I, IV and the X rDNA array. Physical caps and germline continuity are distinct observations, especially at IV 5′ where a cap lies distal to an unresolved gap. The Farm22 output directory is `/lustre/scratch122/tol/teams/blaxter/users/pg17/tol_nematodes/Oscheius_sp_JU75/analysis/nxOscSpeh2.1/ju75_hifi_rescreen.dsGoEq/refine.7W80zD/telomere_clusters.ho2xVx/`. It contains the indexed review FASTA, the standalone X candidate, the alternate tail FASTA, the splice tables, and the supporting read placements.

## I 3′ continuation in the full HiFi BAM

The original coverage BAM has 102 aligned bases at frozen I:8,699,894, the last assembled base. Depth declines from about 500× to about 120× at the internal elimination site and stays above 100× through the tip. Of the reads with MAPQ ≥20 and at least 500 aligned bases distal to the cut, 97 reach the exact tip and 69 carry additional sequence beyond it. These right-hand overhangs are accessible in the full BAM independently of the telomeric read clusters.

The longest continuation is 11,641 bp on `m64230e_230125_115038/13960702/ccs`, with an exact 500 bp splice anchor at MAPQ 60. A second read aligns across 10,593 bp of that continuation at 99.88% identity. Fifty-five other anchored reads reach at least 1 kb into the extension and 15 reach 5 kb. The repeat-rich overhang has many graph matches when aligned on its own. Its aligned proximal sequence establishes the I 3′ placement. The 18 bp of telomeric repeats in the old tip lies inside the extended sequence, and the longest read leaves the distal telomere open.

`JU75.I3.read_extended_candidate.fa.gz` appends the observed 11,641 bp read sequence directly to the review chromosome I with an exact splice. The result is 8,712,803 bp and adds no Ns. The indexed FASTA, read alignment ledger and checksum are in `/lustre/scratch122/tol/teams/blaxter/users/pg17/tol_nematodes/Oscheius_sp_JU75/analysis/nxOscSpeh2.1/ju75_hifi_rescreen.dsGoEq/refine.7W80zD/I3_tip_reads.Fk5BgC/`.

### Comparison with the telomeric end clusters

`C000157`, `C000328` and `C000266` match roughly 3 kb of the new I 3′ sequence at ≥99.2% identity. Their alignments recur within the repeat, with MAPQ 0 for the consensus matches. The most informative member is `C000290` read `E00017617` (`m64230e_230125_115038/122094293/ccs`). It aligns continuously from 820 bp inside the old I 3′ tail through the 11.6 kb extension at 99.34% identity. The last 200 bp before the original tip match exactly and occur once in the frozen reference. The same read continues another 9,012 bp to a 3,949 bp telomeric tract. The tract is 28 bp from the sequenced end. The singleton provides a candidate distal cap on the same molecule. The alignments and `C000290.bridge_evidence.json` are in the `end_clusters.AkPlVr/` subdirectory of the I 3′ run.

## I 5′ distal-cap search

The full HiFi BAM has 39× depth at the first base of frozen chromosome I. After remapping the recruited reads with the other chromosome ends and unplaced sequences present, three high-confidence reads reach the review tip upstream of the frozen start and extend at least 500 bp into the frozen sequence. Two read sequences span the tip and continue outwards. `m64230e_230125_115038/14813941/ccs` supplies 1,312 additional bp with a 500 bp exact splice and MAPQ 51. The other supplies 928 bp overlapping the first. Both end within repeated sequence. Exact, bounded graph walks from the implicated outward ports and the strongest raw-graph placement of the donor read yield zero telomeric paths.

The full telomeric read comparison included 70 kb of I 5′ and the tentative I 3′ sequence as a competing target. `C000157` matches near the extended I 5′ tip and fits the repeat at I 3′. None of its near-tip matches reaches a unique 200 bp anchor in the placed I 5′ sequence. Two `C000251` reads show a shorter telomeric endpoint 4,761 bp inside the frozen I 5′ end. Forty-two recruited reads span that site with at least 500 aligned bases on both sides, including a supported continuation farther out. Two `C000257` reads mark another shorter site 12,591 bp inside the frozen end. The shorter telomeric molecules establish alternate endpoints. Their germline status is unresolved.

The separate, indexed `JU75.terminal_candidates.review.fa.gz` incorporates the 1,312 bp I 5′ extension and the tentative I 3′ `C000290` continuation. It is 57,564,423 bp and carries 11 physical caps under the sequence scanner. I 5′ remains open, I 3′ is tentative, and the IV 5′ cap lies distal to an unresolved gap. The review retains three 200-N gaps. The per-end status, checksums and evidence are beside the FASTA in `/lustre/scratch122/tol/teams/blaxter/users/pg17/tol_nematodes/Oscheius_sp_JU75/analysis/nxOscSpeh2.1/ju75_hifi_rescreen.dsGoEq/refine.7W80zD/I5_tip_reads.Yob2MJ/`.

## Accepted germline caps at both ends of I

**I 3′** is accepted from the `C000290` molecule `E00017617` (`m64230e_230125_115038/122094293/ccs`). The `C000173` family caps the chromosome near the known elimination cut, about 2.4 kb inward of the frozen tip. The accepted read anchors 820 bp into the retained, distal I 3′ tail beyond that cut, including an exact 200 bp stretch with a single occurrence in the frozen reference. The same molecule continues through the 11.6 kb extension, another 9,012 bp, and a 3,949 bp terminal telomeric tract. The curated sequence uses that single contiguous observation, with 28 sequenced bases beyond the annotated tract trimmed. This independently assigns the distal I 3′ port.

**I 5′** is accepted by exclusion under the satellite-port rule. The placed I 5′ tip and the `C000157` cap-facing flank carry the same inward-oriented 224 bp satellite. With the other 11 caps assigned, I 5′ is the remaining open port. The curated sequence prepends the telomeric tract (1,585 bp) and its adjacent 448 bp satellite flank from the single `C000157` molecule `E00013422` (`m64230e_230125_115038/460002/ccs`), in reverse-complement orientation. A 200-N gap between that fragment and the placed tip records the unresolved repeat junction and copy count. The shorter capped `C000251` and `C000257` molecules document alternate endpoints farther inward on the I 5′ arm.

The curated six-chromosome FASTA, `JU75.germline_T2T.curated.fa.gz`, is **57,566,656 bp**. All 12 physical tips pass the canonical telomere scan. The four gaps are 200 Ns each. They are the I 5′ satellite join, the internal I and IV gaps, and the X 5′ rDNA join. The FASTA, index, per-end ledger, manifest, checksums and reproduction script are in `/lustre/scratch122/tol/teams/blaxter/users/pg17/tol_nematodes/Oscheius_sp_JU75/analysis/nxOscSpeh2.1/ju75_hifi_rescreen.dsGoEq/refine.7W80zD/I5_tip_reads.Yob2MJ/curated_T2T.8WCaZe/`.
