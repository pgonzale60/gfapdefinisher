# JU1382 graph and telomere refinement

The JU1382 review FASTA contains six chromosomes totaling **62,601,539 bp**. It adds 353,139 bp relative to frozen `nxOscSpeu1.1`. All **11 frozen N gaps** have a sequence replacement, and the draft contains **zero Ns**. There are **seven supported germline caps**, **one provisional distal cap on IV 5′**, and **four open ends**.

## Inputs

The backbone is the 62,248,400 bp frozen assembly `nxOscSpeu1.1.primary.fa.gz` in the Farm22 freeze directory for `Oscheius_sp._JU1382`. Primary and raw hifiasm unitigs and primary contigs came from `nf_scratch/assemblies/Oscheius_sp_JU1382/results/nxOscSpeu1.pri/hifiasm/`. Original sequence support came from the coordinate-indexed HiFi BAM in `analysis/nxOscSpeu1.1/read_mapping/`, the frozen ONT BAM `analysis/coverage/minimap/ont/nxOscSpeu1.1.bam`, and the 110 assembled telomeric families in `analysis/reads_only/nxOscSpeu1/telomeric_read_fishing/end_clusters/`.

The graph-first `gfapdefinisher finish --nominate-super` run nominated 13 SUPER scaffolds and reported 14 terminal caps. Its 379 output records total 107,640,523 bp. About 40.09 Mb is assigned to the proposed chromosome cores. The curation therefore anchors the three graph exports to the six frozen chromosomes. It checks proposed end paths using original read alignments against the whole frozen assembly, including competing placements. The output adds graph sequence where the required flanks and path are supported. In places where the graph stops, one continuous observed molecule supplies the sequence. ONT intervals are unpolished.

## Chromosome ends

| End | Status | Change in bp | Basis |
| --- | --- | ---: | --- |
| I 5′ | Extended, open | +25,717 | Primary contig `ptg000026l`. Three HiFi molecules enter the extension and two independently anchored ONT reads cover the whole graph segment. The nearby high-depth cap supported by 102 `C000005` reads lies on a shorter branch. Most apparent farther caps anchor principally on II or III. |
| I 3′ | Extended, open | +27,663 | Primary contig `ptg000009l` and three ONT alignments from I into the distal graph sequence. |
| II 5′ | Germline cap | +17,112 | Primary unitig `utg000323l` followed by an ONT allele at a 119 bp exact overlap. Several reads independently anchor in II. An apparent longer branch shares substantial sequence with III 3′. |
| II 3′ | Extended, open | +33,365 | Primary contig `ptg000015l` supported by ONT molecules anchored to II. An internal tract of 199 exact telomeric repeats is about 5.6 kb inward from the draft tip. Multiple reads continue through that tract into the outer sequence. |
| III 5′ | Germline cap | +38,263 | The two unitig graphs and the primary contig agree on the main path. Two ONT reads continue to a distal tract. The chosen read joins the graph at an exact 49 bp anchor. |
| III 3′ | Germline cap | +20,075 | Low-depth three-unitig path represented in both graphs. Two continuous ONT reads independently span from III across the joins to the distal tract. |
| IV 5′ | Provisional distal cap | +3,284 | The outward end of a double-telomeric HiFi molecule supplies the sequence (`C000099`). Its anchored IV interval reaches the catalogued cut at 5,066 bp on the inward side. The case rests on one molecule with additional repeat matches. The longer apparent ONT cap group maps to V. |
| IV 3′ | Open, capped branch and read-supported continuations | 0 | `C000027` and both unitig graphs support a proximal cap. One HiFi reaches the current tip beyond the tract, and two ONT reads continue directly into a partial rDNA monomer. See the [IV 3′ follow-up](#iv-3-follow-up). |
| V 5′ | Germline cap | +3,618 | Primary and raw unitigs provide the same capped sequence. Two `C000047` HiFi molecules and several ONT reads support the terminal junction. |
| V 3′ | Germline cap | +16,376 | Primary unitig `utg000267l` through a 52 bp exact anchor, then one observed ONT allele. The original read carries about 27 kb of V anchor. Its distal repeat also matches II 5′. |
| X 5′ | Germline cap | +12,706 | Identical capped path in both unitig graphs and two original ONT reads with strong X 5′ anchors. Other parts of these molecules match X 3′ through shared sequence. |
| X 3′ | Germline cap | +8,309 | Identical capped representative in both unitig graphs with HiFi and ONT support. One ONT molecule spans the local junction with a strong X 3′ anchor and has a split alignment farther inward. Longer ONT candidates pass through homologous sequence shared between the X ends. |

## IV 3′ follow-up

The 6 bp period in the repeat explorer describes the canonical `TTAGGC` motif. The draft contains 85 exact copies in a 526 bp tract, followed by 494 bp of other sequence. The primary and raw graphs also contain a shorter capped path. Primary `utg000180l` and raw `utg000453l` are identical, with 116 copies and a 4 bp terminal offset. The `C000027` family has 96 HiFi molecules with a unique inward anchor. A somatic interpretation of this cap is plausible given the read count and the paths that continue through the tract.

HiFi molecule `m64016_201224_175337/47581987/ccs` aligns across the tract and the entire 494 bp after it at MAPQ 60. It also supplies the 46 bp patch in the draft. The graph contains the outermost 234 bp at the IV tip. No separate distal telomere was found on this route.

The targeted follow-up screened 97 HiFi and 326 ONT molecules with an inward IV anchor. Two ONT molecules, `0055d8d1-8670-41f3-8c1a-4f14c938a722` and `a6dded7d-d00b-40a1-b517-f90118d7c1f0`, carry separate long IV anchors and then continue through the telomeric tract into the canonical rDNA sequence in forward orientation. Both start the rDNA match at monomer offset 2,065 (zero-based) and share the exact 20 bp junction marker `GGCTTAGGACTACTATTATT`. The alignments span about 3.2–3.4 kb of the monomer. The recruited HiFi reads yielded zero matches to this IV-to-rDNA junction. The ONT molecules end within rDNA and give no distal telomere. Their rDNA route differs from the 494 bp continuation preserved in the draft.

This assigns the proximal cap and supports an additional IV-to-rDNA junction. The distal germline telomere is still open. Follow-up evidence is in `{HiFi,ONT}.IV3.junction_focus.tsv` in the Farm22 run directory.

## Frozen gaps

Coordinates below are zero-based, half-open intervals in the frozen assembly. The sequence-source manifest gives precise replacement intervals, including small flanking trims where needed.

| Frozen gap | Coordinates | Replacement bp | Source |
| --- | --- | ---: | --- |
| I gap 1 | 1,759,377–1,759,577 | 4,206 | HiFi. Several reads span different alleles. |
| I gap 2 | 3,014,127–3,014,327 | 6,518 | ONT. Nineteen molecules bridge the junction in the targeted panel. |
| II gap 1 | 7,371,899–7,372,099 | 1,768 | Primary contig. One HiFi and five ONT molecules cross the tested bridge. |
| II gap 2 | 9,091,167–9,091,367 | 14,361 | Same allele across all three graph exports. Sixteen ONT alignments cross the complete tested rDNA array and at least 2 kb into each non-rDNA flank. |
| IV gap 1 | 1,816,992–1,817,192 | 90,479 | Same allele across all three graph exports. HiFi and ONT alignments span each graph-to-reference boundary. The complete inserted sequence comes from the assembled graph. |
| IV gap 2 | 2,573,830–2,574,030 | 9,575 | Original ONT read `db7bc1df-9fa1-4f2c-9a0a-27139a7e46e1`, uniquely anchored across IV in the whole-reference alignment. The two exact splice flanks are 26 and 30 bp. This is a single unpolished, relatively low-identity allele. |
| IV gap 3 | 6,701,396–6,701,596 | 299 | Same allele across all three graph exports, with 55 exact HiFi sequences. |
| IV gap 4 | 9,677,826–9,677,872 | 46 | Exact HiFi allele. |
| V gap 1 | 30,479–30,679 | 3,383 | Both unitig graphs, supported by HiFi and ONT alignments. |
| V gap 2 | 10,866,928–10,867,128 | 17,794 | ONT. Thirty-one molecules bridge the junction, with modest length variation. |
| X gap 1 | 7,673,179–7,673,379 | 3,063 | Exact HiFi allele, with 11 ONT bridges in the targeted panel. |

The current rDNA array is **internal to chromosome II**. In the freeze, the canonical 7,374 bp monomer and most of a second copy align in the same orientation beginning at II:9,076,780, immediately before gap 2. The graph contributes sequence across that gap. The spanning reads support a continuous internal path and show copy-number variation. The representative sequence comes from the graph. Its repeat count represents one assembly allele. The five archived families with rDNA overlap are singletons and carry no secure distal cap assignment. Four align in the canonical inward orientation (`C000042`, `C000080`, `C000091`, `C000105`). `C000107` aligns in reverse.

## Repeat network

Pairwise end alignments show a shared 13–14 kb subtelomeric block among I 5′, I 3′ and II 3′ at about 94–95% identity. It sits roughly 12 kb inward on I 5′, begins at the current I 3′ draft tip, and sits roughly 16 kb inward on II 3′. The I 3′ draft contains a second homologous copy immediately inward of the first. The largest shared block between III 3′ and II 5′ has about 21.4 kb of matching bases. I 5′, I 3′, II 3′, II 5′ and III 3′ form a network of long homologous blocks. X 5′ and X 3′ share a separate roughly 12.9 kb match. These blocks produce competing end placements and leave several ways to traverse the repeats. The tested reads connect sequence along individual chromosomes. There are zero verified bridges between two distinct chromosome ends. The available repeat-annotation files were empty, so the explorer uses an explicit local tandem screen and reports long sequence similarity separately.

[Interactive chromosome-end explorer](JU1382_satellite_graph/ju1382-satellite-explorer.html) · [GraphML](JU1382_satellite_graph/ju1382-satellite-network.graphml) · [Repeat and read evidence](JU1382_satellite_graph/ju1382-satellite-data.json)

## Outputs

- [Indexed FASTA](../data/benchmark/JU1382/review/JU1382.review.fa.gz) with [FAI](../data/benchmark/JU1382/review/JU1382.review.fa.gz.fai) and [GZI](../data/benchmark/JU1382/review/JU1382.review.fa.gz.gzi).
- [End-status table](../data/benchmark/JU1382/review/JU1382.review.ends.tsv), [sequence-source manifest](../data/benchmark/JU1382/review/JU1382.review.json), and [checksum manifest](../data/benchmark/JU1382/review/JU1382.review.sha256).
- [Best placement of each archived family](JU1382_FAMILY_ASSIGNMENTS.tsv) and [ranked alternatives](../data/benchmark/JU1382/review/JU1382.family_rankings.tsv). Tentative repeat matches retain their lower confidence.
- Farm22 run directory `/lustre/scratch122/tol/teams/blaxter/users/pg17/tol_nematodes/Oscheius_sp_JU1382/analysis/nxOscSpeu1.1/ju1382_graph_refine.XSfJtF/` holds the alignments, original-read subsets, build policy, scripts, and verified assembly.
