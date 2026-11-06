# Biological Classification Taxonomy & PDE Junction Motif Specification

This document defines the biological decision engine specifications for unitig port classification and Programmed DNA Elimination (PDE) junction detection in `gfapdefinisher`.
All classifiers and decision engines defined here strictly adhere to and cite the Biological Truth Rules ([RULES.md](file:///Users/pgonzale60/Documents/research/gfapdefinisher/docs/biology/RULES.md)):
`BIO-001`, `BIO-002`, `BIO-003`, `BIO-004`, `BIO-005`, and `BIO-006`.

---

## 1. The 3-Axis Classification Taxonomy

In a bidirected string graph, every unitig end is modeled as an oriented `SegmentPort` ($u^+$ for 3' suffix, $u^-$ for 5' prefix). Rather than conflating taxonomic origin, local connectivity, and sequence composition into arbitrary ad-hoc labels, `gfapdefinisher` evaluates every port along three orthogonal biological axes:

```
SegmentPort
  │
  ├── 1. Provenance:      nuclear | organelle | contaminant | unknown
  ├── 2. Topology:        tip (0) | collinear (1) | fork (2) | complex_branch (>2)
  └── 3. TerminalFeature: telomeric_cap | rdna_cluster | none
```

### Axis 1: Provenance (Taxonomic & Genomic Origin)
Assigns taxonomic and subcellular compartment provenance based on deep-learning classifications and single-copy orthology:

| Class | Biological Basis & Decision Logic | Rationale & Guardrails |
| :--- | :--- | :--- |
| `nuclear` | `evidence.usco_hit_count > 0` OR `tiara_class == "eukarya"` | Compleasm universal single-copy orthologs (USCOs from `chromadorea_odb12`) or Tiara eukaryotic models confirm genuine nuclear chromosomal sequence. Even if Tiara is uncalled due to sequence length, USCO hit density takes precedence. |
| `organelle` | `tiara_class` in `{"organelle", "mitochondrion", "plastid"}` | High-depth circular mitochondrial or plastid DNA. Must not be incorporated into nuclear telomere-to-telomere scaffolds. |
| `contaminant` | `tiara_class` in `{"bacteria", "archaea", "contaminant"}` | Bacterial or viral contamination originating from host gut microbiome, rearing media, or sequencing reagents. Excluded from eukaryotic walks unless overruled by nuclear USCOs. |
| `unknown` | Unclassified, unassigned, or ambiguous origin | Short fragments or novel sequences lacking sufficient USCO or Tiara taxonomic confidence. |

### Axis 2: Topology (Graph Port Outgoing Degree)
Measures the local connectivity of an oriented port in the bidirected dovetail overlap graph:

| Class | Outgoing Degree | Biological Meaning in String Graph |
| :--- | :--- | :--- |
| `tip` | $\text{deg} = 0$ | Physical dead end in the string graph. Eligible to be an authentic chromosome terminus if capped by a canonical telomeric repeat array (`BIO-001`). |
| `collinear` | $\text{deg} = 1$ | Unambiguous internal continuation along an unbranched chromosomal segment. |
| `fork` | $\text{deg} = 2$ | Simple bifurcation. Characteristic signature of Programmed DNA Elimination diminution boundaries, uncollapsed heterozygous bubbles, or simple repeats. |
| `complex_branch` | $\text{deg} > 2$ | Hyper-branched tangle. Typical of collapsed multivalent repeats, rDNA arrays (`BIO-001`), or satellite clusters. |

### Axis 3: TerminalFeature (Physical Sequence Composition at Terminus)
Classifies the terminal sequence content within the subtelomeric scan window:

| Class | Biological Basis & Decision Logic | Rule Reference |
| :--- | :--- | :--- |
| `telomeric_cap` | `evidence.has_terminal_telomere == True` | Physical end carries a dense, contiguous canonical nematode telomeric repeat array evaluated against `TelomerePolicy` (`min_terminal_repeat_count`, `max_terminal_offset_bp`, `max_internal_gap_bp`) as declared in `src/gfapdefinisher/policy/` as the single source of truth. Authenticates a biological chromosome end. | `BIO-001` |
| `rdna_cluster` | `evidence.has_terminal_rdna == True` | Port enters a 18S-5.8S-28S ribosomal RNA repeat array (Barrnap hit density). An rDNA cluster is an internal repeat tangle, **never** an authentic chromosome end. | `BIO-001` |
| `none` | Neither terminal telomere nor terminal rDNA | Internal continuation, uncapped broken end, or unresolved subtelomeric segment. | `BIO-001` |

---

## 2. Programmed DNA Elimination (PDE) Junction Motif Specification

### Biological Background & Problem Statement
In nematodes undergoing Programmed DNA Elimination (PDE) during somatic differentiation, specific chromosomal regions are excised and eliminated from somatic cell lineages, while being strictly retained in the germline. In sequencing libraries prepared from mixed-stage, whole-worm tissue (e.g. unstaged cultures), DNA from somatic cells ($\sim 95\%$) massively outnumbers DNA from germ cells ($\sim 5\%$).

Consequently, in the uncollapsed Hifiasm string graph (`asm.p_utg.gfa`), somatic diminution boundaries manifest as a characteristic stoichiometric depth drop rather than an assembly error (`BIO-002`, `BIO-003`). Isolated port evaluation cannot differentiate a diminution boundary from a heterozygous bubble or repeat fork. However, topological analysis reveals that genuine PDE junctions follow a universal reachability and stoichiometry motif:

```
                                  ┌───> [ Somatic Cap ] (<= policy.pde.max_somatic_cap_hops,
                                  │                     depth ratio > max_germline_depth_ratio,
[ Somatic Stem ] ─────────────────┤                     telomeric cap)
(e.g. ~64x in reference run)      │
                                  └───> [ Germline Branch ] (depth ratio within
                                        [min_germline_depth_ratio, max_germline_depth_ratio]
                                        relative to core somatic coverage C_core;
                                        e.g. ~4-11x in reference run; continues into GRS; BIO-002)
```
*(Note: Depths ~64x and ~4-11x are concrete empirical values from the DF5120 reference sequencing run; the generalized invariant is the depth contrast ratio declared in policy).*

### Precise Motif Definition
A segment port $p$ is called as a `PDE Junction` if and only if all of the following conditions are satisfied:

1. **High-Depth Somatic Stem**:
   The segment associated with port $p$ has a robust somatic read depth ($rd_{\text{stem}} > 0$).
2. **Bifurcation at Segment Port**:
   Port $p$ has outgoing degree $\text{deg}(p) \ge 2$ (fork or complex branch).
3. **Somatic Branch (Telomeric Cap Reachability)**:
   At least one outgoing branch from $p$ has high read depth ($\text{ratio} = \frac{rd_{\text{branch}}}{rd_{\text{stem}}} > \text{max\_germline\_depth\_ratio}$, declared in `policy.stoichiometry`) AND reaches an authentic terminal telomeric cap within $\le \text{policy.pde.max\_somatic\_cap\_hops}$ traversal hops:
   - **0-Hop (Direct Cap)**: The branch target port itself has outgoing degree 0 (dead-end tip) and carries an authentic canonical telomeric repeat array (`has_terminal_telomere == True`).
   - **Multi-Hop (Linked Cap)**: The branch target port extends through a short linker, adapter, or bubble unitig ($\le \text{max\_somatic\_cap\_hops}$) to an authentic dead-end telomeric cap tip with somatic-level depth.
4. **Germline Branch (Stoichiometric GRS Continuation)**:
   Concurrently, at least one outgoing branch from $p$ exhibits a stoichiometric depth ratio relative to core somatic coverage ($C_{\text{core}}$) in the calibrated range:
   $$\text{policy.stoichiometry.min\_germline\_depth\_ratio} \le \frac{rd_{\text{germline}}}{C_{\text{core}}} \le \text{policy.stoichiometry.max\_germline\_depth\_ratio}$$
   governed strictly by `policy.stoichiometry` (`BIO-002`). Candidate germline branches are ranked by biological criteria (proximity to expected GRS stoichiometry $( \text{min} + \text{max} ) / 2$, segment length, and branch depth support), using unitig port name strictly as a deterministic tie-breaker. The germline continuation branch must be distinct from the somatic terminus branch.
5. **Complex GRS Retention (`BIO-004`, `BIO-005`)**:
   A candidate germline branch is accepted based on depth stoichiometry contrast regardless of whether it continues as an unbranched segment or leads into a complex repeat, satellite array, or graph tangle. Retaining the initial GRS entry from a confident PDE junction ensures partial GRS sequences and satellite arrays are preserved (`BIO-004`) rather than falsely discarding the diminution call.
6. **Deterministic Label Assignment**:
   When detected, the junction call deterministically outputs:
   - `somatic_port`: The outgoing branch port heading into the somatic telomeric terminus.
   - `germline_port`: The outgoing branch port continuing into eliminated germline-restricted sequence (GRS).
   - `depth_ratio`: The empirical germline-to-core depth ratio ($D / C_{\text{core}}$).
   - `stem_depth_ratio`: The empirical germline-to-stem depth ratio ($D / rd_{\text{stem}}$).
   - `core_depth_ratio`: The empirical germline-to-core depth ratio ($D / C_{\text{core}}$).

---

## 3. Negative Cases & Falsification Guardrails

To prevent false-positive diminution calls, the classifier explicitly rejects:

1. **Symmetric Heterozygous Bubbles**:
   Bifurcations where both outgoing branches have comparable read depths ($\sim 0.50 \times C_{\text{core}}$) and neither terminates in an authentic telomeric cap.
2. **Spurious Low-Depth Artifacts / Chimeras**:
   Branches with depth ratio below $\text{policy.stoichiometry.min\_germline\_depth\_ratio}$ (declared in `policy.stoichiometry`), representing sequencing noise, coverage dropout, or chimeric reads rather than true germline retention.
3. **Internal Collapsed Repeats**:
   Forks where branches exceed $\text{policy.stoichiometry.max\_germline\_depth\_ratio}$ without terminating in a telomeric cap (e.g. internal transposon or satellite duplications).
4. **Uncapped Somatic Termini**:
   Forks where a putative somatic branch lacks verified telomeric repeats (`BIO-001`).

---

## 4. Benchmark Validation & Generalization Guard

In accordance with `BIO-006`, the PDE junction motif detector has been calibrated and validated against reference diminution boundaries in benchmark dataset `DF5120` (where length-weighted median core coverage $C_{\text{core}} = 63.0\times$):

| Benchmark Junction | Stem Port | Somatic Branch Port | Somatic Depth (Ratio) | Germline Branch Port | GRS Depth (Ratio to C_core) | Status |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| `utg000001l:3p` | `utg000001l+` | `utg000142l-` | $46.0\times$ ($0.7302$) | `utg000101l+` | $6.0\times$ ($0.0952$) | `VALIDATED` |
| `utg000001l:5p` | `utg000001l-` | `utg000076l+` | $27.0\times$ ($0.4286$) | `utg000039l+` | $7.0\times$ ($0.1111$) | `VALIDATED` |
| `utg000002l:3p` | `utg000002l+` | `utg000194l+` | $47.0\times$ ($0.7344$) | `utg000115l+` | $8.0\times$ ($0.1270$) | `VALIDATED` |
| `utg000002l:5p` | `utg000002l-` | `utg000093l-` | $48.0\times$ ($0.7500$) | `utg000103l+` | $11.0\times$ ($0.1746$) | `VALIDATED` |
| `utg000009l:5p` | `utg000009l-` | `utg000211l-` | $27.0\times$ ($0.4355$) | `utg000069l+` | $6.0\times$ ($0.0952$) | `VALIDATED` |
| `utg000013l:3p` | `utg000013l+` | `utg000167l+` | $60.0\times$ ($0.9524$) | `utg000033l-` | $8.0\times$ ($0.1270$) | `VALIDATED` |
| `utg000027l:5p` | `utg000027l-` | `utg000154l-` | $60.0\times$ ($0.9231$) | `utg000189l-` | $8.0\times$ ($0.1270$) | `VALIDATED` |
| `utg000016l:5p` | `utg000016l-` | `utg000088l+` | $45.0\times$ ($0.7143$) | `utg000096l+` | $7.0\times$ ($0.1111$) | `VALIDATED` |
| `utg000019l:3p` | `utg000019l+` | `utg000200l-` | $39.0\times$ ($0.5909$) | `utg000071l+` | $7.0\times$ ($0.1111$) | `VALIDATED` |
| `utg000019l:5p` | `utg000019l-` | `utg000129l-` | $37.0\times$ ($0.5606$) | `utg000145l+` | $4.0\times$ ($0.0635$) | `VALIDATED` |
| `utg000021l:5p` | `utg000021l-` | `utg000051l+` | $51.0\times$ ($0.7969$) | `utg000071l-` | $7.0\times$ ($0.1111$) | `VALIDATED` |
| `utg000027l:3p` | `utg000027l+` | `utg000209l+` | $66.0\times$ ($1.0154$) | `utg000145l-` | $4.0\times$ ($0.0635$) | `VALIDATED` |

Per `BIO-006`, this decision logic uses pure topological connectivity and depth ratios without hardcoded unitig identifiers or species-specific parameters, ensuring unconditional transferability across held-out outcrossers (`JU1382`, `PS1017`, `PS2068`).
