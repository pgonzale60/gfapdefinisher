# DF5120 germline terminal evidence

## Reproduce the milestone

The accepted DF5120 assembly contains six germline telomere-to-telomere scaffolds, 12 assigned distal caps and 48,227,278 bp. It has ten 200 N joins across two regions. Nine are in chromosome II rDNA and one is at III 5′. A 200 N run marks an uncertain connection with no inferred physical length.

The [DF5120 recipe](../benchmarks/df5120_t2t.json) records reviewed graph paths and read splice coordinates. `finish` builds the graph-based scaffolds. `curate` uses those and the primary and raw GFAs to replay the selected paths and extract the two additional ONT bridges directly from the reads. The prior DF5120 reference supplies chromosome names and orientation. The recipe runs without the historical checkpoint or the accepted FASTA.

From the repository root with the DF5120 benchmark files staged, run:

```bash
data=data/benchmark/DF5120
out=data/benchmark/DF5120/milestone-rebuild

uv run gfapdefinisher finish \
  --gfa "$data/asm.p_utg.gfa.gz" \
  --tiara "$data/tiara_unitigs.txt" \
  --barrnap "$data/barrnap_rdna.gff" \
  --ont-paf "$data/ont_to_utgs.paf" \
  --ont-reads "$data/ont_spanning_reads.fa.gz" \
  --nominate-super --out-dir "$out/finish"

uv run gfapdefinisher curate \
  --recipe benchmarks/df5120_t2t.json \
  --scaffolds "$out/finish/scaffolds.fasta" \
  --gfa "$data/asm.p_utg.gfa.gz" \
  --raw-gfa "$data/asm.r_utg.gfa.gz" \
  --ont-reads "$data/ont_spanning_reads.fa.gz" \
  --out "$out/DF5120.fa"
```

The expected FASTA SHA-256 is `f6b95df1eebde82009db7026308a4447e5103c3ba0a60dcced8c65d711c7482d`. The command checks every declared exact overlap and both splice anchors for each read bridge, confirms the assigned cap sequence and terminal tract at each end, and checks the gap counts and output checksum before writing the FASTA.

Two internal 200 N gaps in the frozen reference are covered in the automatic `finish` output and in the milestone. I contains the 5,472 bp graph interval, and IV contains the 3,271 bp graph interval. In each case the primary and raw GFA intervals are identical, and the full interval plus its two flanks was verified in both outputs. `finish` builds those sequences from the graph without comparing against the reference FASTA. The frozen internal V gap falls in a shorter graph representation, with no verified replacement for the missing reference interval. The two regions still marked with Ns in the milestone have no verified graph-only fill. [Frozen-reference gap comparison](DF5120_GRAPH_REFINEMENT.md#published-assembly)

Use the existing chromosome assembly to orient and place graph paths. A recipe may also use existing sequence with `--reference`. Where the reference and the reviewed GFA path disagree, select the GFA path in the recipe. The IV and V edits, for example, extend past the somatic telomeric tract to the distal germline cap. Reads from HiFi or ONT can supply the sequence between anchored flanks. One continuous supported read sequence is sufficient for a satellite join and leaves zero Ns. The selected DF5120 bridges use ONT and retain the chosen read sequence without polishing.

When the repeat and available ports support a unique placement, record the distal cap across an unknown join if sequence spanning the connection is unavailable. The DF5120 III 5′ placement remains based on curated port inference. Keep the available rDNA copies in a consistent transcriptional orientation between the assigned flanks and retain 200 N markers at uncertain joins. The milestone replays these reviewed decisions. Discovery and adjudication of every port in another isolate are separate from replay.

## Historical evidence and compilation

Earlier curation identified a distal telomere for all 12 chromosome ends when the evidence from the different runs is considered together. IV 3′ and V 3′ have ONT evidence for placing their distal caps across variable satellite sequence. III 5′ and the chromosome II rDNA end rely on scaffold placement across repeats. The historical FASTAs checked in this audit each contain a subset of those placements.

## Compiled FASTA

`data/benchmark/DF5120/nxOscSper1.DF5120.distal_caps.fa` combines the 12 assigned caps in six chromosome records, oriented as in `nxOscSper1.1`. It contains 48,262,687 bp and has SHA-256 `e914c17f39fc11d9c4de8a4d47f1cf816736554fbc7939001e23fa0094ea3f49`. Both terminal 2 kb sequences on each chromosome match the corresponding primary-graph cap in the orientation listed below.

I, III and V come from the staged checkpoint. II and X are the reverse complements of `SUPER_2` and `SUPER_1` from the `796d9f1` run with staged ONT inputs. IV takes the first 8,356,957 bp of `SUPER_3` and the checkpoint suffix beginning after checkpoint IV base 8,279,657. The shared core sequence through the splice matches exactly. This keeps the ONT-supported IV 5′ end and the historical gapped IV 3′ cap placement. The raw IV 3′ extension remains a candidate for a future sequence-level edit.

The FASTA retains nine 100 N gaps in II rDNA, an 11,000 N gap at III 5′, an 8,500 N gap at IV 3′ and a 5,000 N gap at V 3′. The III 5′ placement remains provisional. At MAPQ ≥ 20 on the expected chromosome, alignments cover 26,493 of 69,446 curated IV 3′ GRS bases and all 18,754 curated V 3′ GRS bases.

### Read-bridged FASTA

`data/benchmark/DF5120/nxOscSper1.DF5120.read_bridged.fa` uses one continuous ONT sequence at each of the IV 3′ and V 3′ satellite joins. The cuts have exact 24 bp read-to-assembly matches centered on each splice point. There was no consensus or polishing, so the inserted sequence represents the selected read. Coordinates below are zero-based, half-open. The read intervals are in chromosome orientation after reverse complementing the named reads.

| End | Replaced interval in compiled FASTA | ONT read | Inserted interval |
| --- | --- | --- | --- |
| IV 3′ | `8,356,948–8,377,510` | `e87d54ff-c85d-4639-96a0-595eb018f886` | `8,517–11,235` (2,718 bp) |
| V 3′ | `7,027,827–7,036,371` | `ac008920-a981-414e-bc11-3a2dcb6ea4c8` | `33,628–34,507` (879 bp) |

Both source reads align continuously across their new joins at MAPQ 60. All 12 terminal 2 kb sequences match the compiled FASTA exactly. The read-bridged FASTA contains 48,227,278 bp with SHA-256 `f6b95df1eebde82009db7026308a4447e5103c3ba0a60dcced8c65d711c7482d`. Its nine rDNA gaps on II and one satellite gap on III are each exactly 200 Ns. The Ns mark uncertainty in the connections, independently of gap length. The ten full HiFi reads used for the earlier III consensus align best within individual unitigs when remapped against both whole graphs. Their repeat-only split alignments to the old flank references are insufficient to place the patch.

The user-provided frozen ONT alignment is `freeze/2022-02-17/data/Oscheius/Oscheius_sp._DF5120/analysis/coverage/minimap/ont/nxOscSper1.1.bam` on Farm22. All 374 primary read records intersecting `III:1-80000` were remapped on Farm22 against both graphs. None gave separate, ordered alignments in the expected orientation from the assigned cap into the core-side unitig or its raw-graph GRS candidate. The window screen covers reads mapped to that end of the frozen assembly. A read with a different original placement could fall outside it.

### Chromosome III satellite k-mers

In the terminal 2 kb of the chromosome-side unitig `utg000103l`, exact 15-mer matches recur most often at offsets of 132 and 660 bp in the tested ranges of 120–145 and 630–690 bp. The 660 bp spacing is consistent with five 132 bp copies. This measures periodicity in the assembled repeat. Of the 222 distinct canonical 31-mers occurring at least three times in that window, 222 occur near the telomere on `utg000150l` and zero occur in its gap-facing 2 kb. The repeat-facing ends of raw `utg000184l` and `utg000476l` both contain all 222. The gap-facing 5 kb of the cap occurs exactly in reverse orientation inside primary `utg000083l` and raw `utg000184l`. The raw graph links `utg000184l-` to `utg000476l+`. Together these observations identify a candidate intervening path and explain why a match to the 132 bp family by itself cannot place the join.

The Farm22 LSF read screen (array `50451`) uses exact 21-mers unique across the complete primary graph. It found 45 such markers on the cap, 1,202 on `utg000083l` and 3,204 within 40 kb of the outward end of `utg000103l`. A candidate needs separate, ordered chains of at least three non-overlapping markers on both pieces. This screen is sensitive to sequencing errors in the small set of exact cap markers. Farm22 run directory is `/lustre/scratch122/tol/teams/blaxter/users/pg17/nf_scratch/assemblies/Oscheius_sp_DF5120/analysis/iii5_satellite_kmers.NNV1Va`.

Both LSF jobs finished successfully. The screen found zero candidates among 503,743 HiFi reads and zero among 922,991 ONT reads. The III 5′ join remains represented by a 200 N gap.

### Estimating satellite abundance from reads

The separate Farm22 abundance screen (array `53662`) counted 202 recurring canonical 21-mers from the III core-side repeat and independently counted graph-unique 21-mers from the assembled cap, intermediate unitig and germline flank. Exact repeat matches starting one base apart were tallied as adjacent hits. For sizing, overlapping 21-mer intervals on a read were merged so adjacent hits contribute only the bases that they add. The 202 markers occur only on III in the frozen assembly, covering 36,074 bp there and 19,725 bp in the read-bridged FASTA. In the primary graph they are concentrated in `utg000083l` (62,053 bp), `utg000103l` (18,732 bp) and `utg000150l` (993 bp). Those unitig footprints include unplaced sequence and cannot be added to an assembled chromosome without determining their order and overlap.

Across 503,743 HiFi reads, 107 carried at least 20 repeat hits. They contained 729,823 occurrences, including 724,621 whose start is one base after the previous repeat hit. Their nonoverlapping repeat footprint was 825,143 sequenced bases. Dividing by 11×, the median exact-marker depth of the assembled germline flank, gives approximately **75 kb**, or **568 132 bp monomer equivalents**. The cap markers give 7× and an alternative normalization of approximately **118 kb**, or **893 monomer equivalents**. The retained chromosome III core gives 74× and is a control for depth contrast. The two germline depths show sensitivity to the choice of denominator. This estimates the abundance of the recognized family, including assembled copies. It does not measure the length of the missing junction.

Seven HiFi reads carried both repeat markers and cap markers, and 38 carried repeat markers and germline-flank markers. All seven cap-associated reads had repeat footprints of about 1 kb, consistent with the array already represented near that cap's telomere. In all seven cap-associated reads and 28 germline-flank-associated reads, at least one flank marker lay outside the read's outermost repeat-marker hits. None of the HiFi reads carried both sets of flank markers.

Across 922,991 ONT reads, 160 met the same 20-hit threshold. They carried 420,752 repeat hits, including 410,175 adjacent starts, and 625,419 bp of merged repeat footprint. The 9× germline-flank depth gives **69.5 kb**, or **526 monomer equivalents**. The 8× cap depth gives **78.2 kb**, or **592 equivalents**. Two ONT reads met the marker-count threshold for both the cap and germline flank, but their cap markers fall inside the overall repeat-marked interval. These are repeat-variant matches without independent outside anchors on both sides. Taken together, the germline-flank normalizations give a working recognized-family abundance of **70–75 kb**, about **530–570 monomer equivalents**. Both LSF elements finished. Farm22 run directory is `/lustre/scratch122/tol/teams/blaxter/users/pg17/nf_scratch/assemblies/Oscheius_sp_DF5120/analysis/iii5_satellite_abundance.IlC9q8`.

## Which outputs the records describe

- The staged `nxOscSper1.germline_T2T_validated.fa` is byte-for-byte identical to the Farm22 `analysis/hifiasm_t2t_validation` output (SHA-256 `401617743e5bb03622a0ee79d084a86c0c03a2f36ab56b6852f75b4a7eea00fe`). Its `junctions.tsv` and `unbridged_gaps.bed` place telomeric unitigs across three satellite gaps. The FASTA stops inside rDNA at II 3′ and has a 6.81 Mb X that terminates with a somatic tract at one end. The archived command for this exact FASTA used the primary GFA and both unitig FASTAs, with no spanning-read patch input. The work leading to that FASTA included separate read analyses. [Original toolkit #71](https://github.com/pgonzale60/PDE_origin/pull/71), [checkpoint changes #74](https://github.com/pgonzale60/PDE_origin/pull/74)
- The merge snapshots for gfapdefinisher #63 (`47c098a`) and #64 (`42e1e97`) were rerun with the staged primary GFA, Tiara, TIDK, Barrnap, Compleasm, ONT PAF and ONT FASTA. Both reruns and the current `796d9f1` output place the II rDNA cap and build the 9.26 Mb X. All three outputs end on the core-unitig somatic cap at IV 3′ and V 3′. The regression test checks motif abundance, length and selected components. It explicitly accepts the `utg000009l` core-cap end state. [gfapdefinisher #63](https://github.com/pgonzale60/gfapdefinisher/pull/63), [gfapdefinisher #64](https://github.com/pgonzale60/gfapdefinisher/pull/64), [current regression test](../tests/test_cli.py)
- The separate PDE_origin curation dossier describes an `nxOscSper1.2.rc2.primary.fa` with all 12 ends and one internal V gap. The named FASTA was absent from the local files and the Farm22 assembly and staging directories searched. The dossier has two different chromosome II lengths in its assembly summary and AGP excerpt, so its all-end claim could not be checked against a specific FASTA. [Curation dossier](https://github.com/pgonzale60/PDE_origin/blob/21a194b/docs/curation_dossiers/nxOscSper1_DF5120_curation_dossier.md)

## Evidence at each end

End labels use published `nxOscSper1.1` orientation. In current `SUPER` output, II, V and X are reversed. Unitig IDs in the table refer to the primary graph unless marked raw.

| End | Distal cap | Core-to-cap evidence and remaining uncertainty |
| --- | --- | --- |
| I 5′ | `utg000039l` | Exact 15,643 bp overlap with the core in the Farm22 junction ledger. |
| I 3′ | `utg000101l` | Exact 14,857 bp overlap with the core in the same ledger. |
| II 5′ | `utg000205l` | Graph path via germline `utg000033l`. The core join is exact. The cap join has lower sequence identity in the Farm22 ledger. |
| II 3′ | `utg000053l` | The core enters rDNA through a graph edge to `utg000020l`. The distal cap is a separate telomeric rDNA unitig selected by orientation and rDNA content. #64 placed oriented copies and marked uncertain joins with AGP gaps. The internal array and copy count remain unresolved. [Issue #26](https://github.com/pgonzale60/gfapdefinisher/issues/26), [PR #64](https://github.com/pgonzale60/gfapdefinisher/pull/64) |
| III 5′ | `utg000150l` | The core-to-GRS join to `utg000103l` is exact. The cap is placed across a repeat gap. An earlier run proposed a 3,961 bp HiFi consensus from ten reads. Rechecking the staged flank PAF, none of those ten had MAPQ ≥ 20 on the repeat-side flank. The subsequent recipe restored the gap. This is the least secure distal assignment. [PDE_origin #72](https://github.com/pgonzale60/PDE_origin/pull/72), [PDE_origin #74](https://github.com/pgonzale60/PDE_origin/pull/74), [gfapdefinisher #65](https://github.com/pgonzale60/gfapdefinisher/pull/65) |
| III 3′ | `utg000115l` | Exact 4,706 bp overlap with the core in the Farm22 ledger. |
| IV 5′ | `utg000217l` | gfapdefinisher puts the cap beyond `utg000069l` using two ONT reads and emits 1,187 bp from one unpolished read. The rerun recovers the full curated GRS by alignment span at this end. The old staged checkpoint had rejected the graph-only 5′ join. [gfapdefinisher #63](https://github.com/pgonzale60/gfapdefinisher/pull/63) |
| IV 3′ | `utg000056l` | The staged checkpoint uses an 8,500 N gap after trimming the somatic cap. Two archived ONT reads independently align in the expected order and orientation from the core flank to the distal unitig at MAPQ 60 on both sides. The raw graph adds an anchored GRS continuation. Repeat copy number and the exact sequence join remain uncertain. |
| V 5′ | `utg000139l` | Exact 3,852 bp overlap with the core in the Farm22 ledger. |
| V 3′ | `utg000126l` | The staged checkpoint trims the somatic cap, extends through raw `utg000309l` carrying the entire curated GRS, and places the distal cap across 5,000 Ns. One archived ONT read maps at MAPQ 60 to both the raw continuation and distal raw unitig `utg000274l` in the expected order and orientation. The alignments overlap in the repeat, leaving the precise join and gap size uncertain. [PDE_origin #74](https://github.com/pgonzale60/PDE_origin/pull/74) |
| X 5′ | `utg000096l` | The longer primary-graph walk restores the published arm, including the curated terminal and internal GRS. It reaches the distal cap. Several graph overlaps need sequence-level review. [Graph refinement audit](DF5120_GRAPH_REFINEMENT.md) |
| X 3′ | `utg000189l` | Core-to-cap graph walk, supported by a 3,900 bp exact junction in the Farm22 ledger. |

## Read checks on the two former somatic caps

The Farm22 `ont_gap_assessment` archived 27 ONT candidates. They were realigned against the whole primary and raw unitig FASTAs with minimap2 `-x map-ont -c` and assessed for orientation, query order, and alignment on both sides of each target gap. The IV evidence is on reads `ad4bdee6-c279-4880-8033-83137b8ac0ae` and `e87d54ff-c85d-4639-96a0-595eb018f886`. Their primary-graph distal alignments span 11.1 and 9.4 kb at MAPQ 60. The V read `ac008920-a981-414e-bc11-3a2dcb6ea4c8` aligns to bases 5,740–43,009 of raw `utg000309l` and bases 50,358–55,434 of raw `utg000274l` in opposite target orientations, as expected for this placement. On the read the alignments overlap by 2,646 bp in the repeat. An old flank set used the somatic stub `utg000153l` for V, so support counted against that target cannot validate the distal V cap.

Remapping all 15,838 reads in the staged ONT subset reproduced the IV and V alignments. A screen requiring collinear alignments of at least 1 kb at MAPQ ≥ 20 on both sides yielded no support for the III 5′ join in that subset.

For a benchmark of chromosome-end placement, the ledger needs to track the expected distal cap and the support for reaching it by each route. The staged validation FASTA and the motif test alone give different answers because they preserve different parts of the curation.
