# PS1017 satellite connections

Open [the explorer](ps1017-satellite-explorer.html) in a browser. The graph loads locally and can be explored offline. The [GraphML network](ps1017-satellite-network.graphml) carries the same relationship classes, array intervals, read IDs, exact k-mer counts, and junction-site attributes for network software such as Cytoscape. [Cytoscape network formats](https://manual.cytoscape.org/en/latest/Supported_Network_File_Formats.html)

## Reading the graph

- Click a **satellite group** for the chromosome ends and fragments that carry it, array spans, periodicity, junction counts, and oriented 21-mer comparisons. The estimated spans include a 300 kb survey at each assembled chromosome end.
- Click a **chromosome end, fragment, or graph path** for its repeat track, recorded junctions, and read evidence. Solid green edges represent observed molecules across a particular join. Gray-blue edges represent local repeat occurrence. Dashed amber edges represent competing placements or proposed graph paths.
- **Junctions** switches between distinct junction sites and all 56 archived read families. The catalog records seven placed distal caps, 26 distinct catalogued somatic sites represented by these families, and one observed cap on an unplaced fragment. C000019 contributes two site records at IV 3′. Candidate placements remain listed with their confidence or supporting context.

The [static figure](ps1017-satellite-figures.png) shows array coverage across the locus set and detailed tracks for II 5′, III 5′, the competing V path, and both ends of X. [Model data](ps1017-satellite-data.json)

## Repeat groups and junctions

The mixed 175–343 bp group is the broadest in this survey. It occurs at both X ends, the revised II 5′ tip, I and V, and several frozen fragments. X 5′ has a predominantly 343 bp array. X 3′ and the new II 5′ tip contain mixed variants. Exact 21-mer comparisons help distinguish these variants within the broader group.

The **23 bp repeat** has two surveyed records. The first is the array now crossed by the gap-free II 5′ splice. The second is the corresponding edge of frozen `II_unloc_3`, whose sequence was incorporated into II. Four original ONT molecules span that join. Both labels are preserved for provenance and refer to the same incorporated path.

The **20 bp repeat** exposes the unresolved III 5′ and V 3′ placement. C000038 supports a telomere on the capped graph fragment `utg000068l`. Its opposite edge contains the short repeat and matches both the III port and the longer V 3′ graph path. The fragment contributes one observed junction, with two possible chromosome placements. Its repeat and cap occur at different positions along that fragment.

“Adjacent” counts junctions within 1 kb of a surveyed array. “On these loci” also includes greater separations such as the junction on the C000038 fragment. These are survey counts for the reviewed assembly and selected evidence, with additional unresolved families listed separately.

## Sources and coordinates

The assembly is `ii5_chain_extension/PS1017.II5_extended.review.fa.gz` in the [Farm22 run directory documented here](../PS1017_GRAPH_REFINEMENT.md). The model combines the frozen satellite-array GFF3, the prior family rankings and the audits of spanning reads and graph placements. The eight repeat groups cluster overlapping family labels where appropriate and include the audited 20 bp and 23 bp arrays.

Coordinates are zero-based and half-open. Chromosome tracks are measured inward from the current tip. Frozen fragment tracks start at the native 5′ end, and graph tracks use the proposed outward edge. Lifted TRF and marker-based spans carry separate evidence labels. Sections of newer graph paths lacking TRF annotation are surveyed with recurrent 11-mer markers. An exact 21-mer comparison samples up to three of the largest intervals for each locus and family, capped at 25 kb per interval.

The native hifiasm graph record names are visible in the explorer. Examples include `utg000068l` for the C000038 cap fragment, `utg000063l` for the farther II path, and `ptg000009l` for the longer V 3′ path. The GraphML edges encode these relationships for exploration. The native GFA retains the hifiasm sequence graph.
