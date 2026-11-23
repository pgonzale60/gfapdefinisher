# DF5120 graph refinement audit

The HiFi graphs support two gap fills in the published assembly and a 2.44 Mb extension to the staged checkpoint's chromosome X. The raw graph also recovers the core connection to germline sequence on IV and V. The three repeat gaps in the checkpoint remain open.

## Inputs and comparison

- Published `nxOscSper1.1`, six nuclear chromosomes from GenBank accession `GCA_932521415.1`. Chromosomes were named I, II, III, IV, V and X for the comparison. [GenBank FASTA](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/932/521/415/GCA_932521415.1_nxOscSper1.1/GCA_932521415.1_nxOscSper1.1_genomic.fna.gz)
- Staged `data/benchmark/DF5120/nxOscSper1.germline_T2T_validated.fa`. This is a separate assembly with different terminal sequence and gap positions.
- Staged `asm.p_utg.gfa.gz` and `asm.r_utg.gfa.gz`, with curated GRS from `nxOscSper1.1_GRS.fasta` used to identify germline sequence. Unitig IDs are specific to each graph.

The Farm22 freeze copy at `/lustre/scratch122/tol/teams/blaxter/users/pg17/freeze/2022-02-17/data/Oscheius/Oscheius_sp._DF5120/assembly/curated/nxOscSper1.1.primary.fa.gz` matches the six GenBank nuclear chromosomes base for base after normalizing FASTA headers, line wrapping and case. Both copies contain 47,746,511 bp and the same three 200 bp N-gaps.

Extract the GFA `S` sequences and align them to each FASTA with minimap2 2.31-r1302 (`-x asm5 -c`). Align up to 50 kb on each side of every N-run of at least 5 bp back to each graph's unitig FASTA. Check strand, query order, CIGAR, GFA links and concordance between the graphs. Align curated GRS against the unitigs independently. Gap coordinates below are zero-based, half-open. End labels are in the respective FASTA's orientation.

## Published assembly

| Gap | Primary and raw graph evidence | Result |
| --- | --- | --- |
| I `3,778,464–3,778,664` | Both 50 kb flanks map exactly at MAPQ 60, in order, to primary `utg000001l+` and raw `utg000074l+`. The intervening sequence is identical across both graphs and the checkpoint. | Replace the 200 Ns with **5,472 bp** of graph sequence for this graph-supported representation. The two graphs use the same underlying HiFi data. |
| IV `6,351,203–6,351,403` | Both 50 kb flanks map at MAPQ 60, in order, to primary `utg000009l+` and raw `utg000024l-`. The intervening sequence is identical across both graphs and the checkpoint. | Replace the 200 Ns with **3,271 bp**, a net gain of **3,071 bp**. |
| V `1,958,669–1,958,869` | Primary `utg000004l-` and raw `utg000076l-` align on both sides. The separately mapped 50 kb flanks overlap by 37 bp in graph coordinates, with different sequences in those 37 bp on the reference. The long alignments absorb 13,742 reference bases in a single deletion across this locus. | Keep the gap. The graph presents a shorter repeat representation and provides no unambiguous insert between the existing flanks. |

The primary-graph insert intervals are `utg000001l:3736757–3742229` for I and `utg000009l:6278270–6281541` for IV, both on the forward strand. The I interval lies in a tandem repeat. The local `nxOscSper1.2` curation dossier describes a different, 10 bp overlap at this gap, so the graph-supported insert is a candidate representation and the repeat copy count needs review before adoption.

The primary graph also has seven end alignments with a telomeric overhang. Approximate overhangs from the alignments are listed below. These quantify additional terminal sequence, which can include telomeric repeats.

| Published end | Primary unitig | Overhang |
| --- | --- | ---: |
| I 5′ | `utg000039l-` | 17,817 bp |
| I 3′ | `utg000101l+` | 11,051 bp |
| III 3′ | `utg000115l+` | 10,801 bp |
| IV 3′ | `utg000056l+` | 3,755 bp |
| V 5′ | `utg000139l-` | 7,658 bp |
| X 5′ | `utg000096l-` | 8,664 bp |
| X 3′ | `utg000189l-` | 2,775 bp |

At published V 3′, raw `utg000309l+` matches the final 28,971 bp exactly and extends another **14,042 bp** without reaching a telomere. It contains all 18,754 bp of the curated terminal GRS and connects to the core through `utg000111l-`. The existing checkpoint already includes the full raw unitig.

## Staged checkpoint

**X 5′** has the largest extension candidate. The complete 6,814,393 bp checkpoint X matches one end of the primary-graph chromosome sequence exactly in reverse orientation. The graph continues **2,442,336 bp** through six additional unitigs to a telomeric tract. It recovers the published X arm, including the terminal GRS and the internal GRS near published position 2.34 Mb. The published assembly already includes this arm.

The linked path out of the checkpoint toward X 5′ is shown below in graph orientation. The 2,442,336 bp estimate trims incoming GFA overlaps. Three of the six overlap sequences differ when compared directly, so the path is an extension candidate and the exact join sequences need review.

```text
utg000019l+ → utg000071l+ → utg000021l+ → utg000081l+ → utg000065l+ → utg000016l- → utg000096l+
```

| Checkpoint gap | Graph evidence | Result |
| --- | --- | --- |
| III 5′, 11,000 Ns | The two flanks map exactly to distinct unitigs in both graphs. Raw `utg000476l` recovers part of curated GRS and aligns near the core-side flank. The exact telomere-side unitig has no outgoing link toward the gap. | Keep the gap. |
| IV 3′, 8,500 Ns | Raw `utg000092l+` links to `utg000160l-`. That unitig matches 11,399 bp immediately before the gap exactly and supplies 15,572 bp beyond the flank. Its match to the distal flank lies in repetitive sequence. The distal telomere-bearing `utg000132l+` has no graph link into this walk. | The raw graph extends from the core into GRS. Keep the gap until the distal join can be placed. |
| V 3′, 5,000 Ns | Raw `utg000111l-` links to `utg000309l+`, which maps exactly up to the gap. The distal telomere-bearing unitig `utg000274l-` is disconnected. | The raw graph confirms the core-to-GRS placement already present in the checkpoint. Keep the gap. |

At checkpoint II 3′, the core links to the first rDNA-bearing unitig in the primary graph. The existing FASTA already reaches into rDNA. The rest of the array remains a gapped scaffold between its identified flanks under the rDNA policy. None of the candidate edits above requires linking floating GRS scaffolds as a separate step. This comparison uses graph sequence and topology alone, and the input FASTAs were kept intact.
