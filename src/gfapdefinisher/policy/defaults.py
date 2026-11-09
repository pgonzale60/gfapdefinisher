"""Default policy declarations with mandatory justifications and provenances."""

from dataclasses import dataclass

from gfapdefinisher.policy.models import PolicyValue


@dataclass(frozen=True)
class TelomerePolicy:
    motif_forward: PolicyValue[str] = PolicyValue(
        value="TTAGGC",
        justification="Canonical nematode telomeric repeat hexamer across Rhabditida.",
        provenance="Nematode telomere literature standard; BIO-001",
    )
    motif_reverse: PolicyValue[str] = PolicyValue(
        value="GCCTAA",
        justification="Reverse complement of canonical nematode telomeric repeat hexamer.",
        provenance="Mathematical derivation; BIO-001",
    )
    min_terminal_repeat_count: PolicyValue[int] = PolicyValue(
        value=30,
        justification="we want as low a number possible to still consider telomeric reads that did not span the whole telomere, but large enough to discrimate internal telomere sequences. It is open for revision upon bigger data observation",
        provenance="@pgonzale60 biological calibration on reference benchmark and nematode telomere architecture (BIO-001)",
        provisional=True,
    )
    terminal_scan_bp: PolicyValue[int] = PolicyValue(
        value=10000,
        justification="Maximum sequence length inspected from segment physical terminus for telomeric arrays; BIO-001",
        provenance="Subtelomeric terminal window standard",
        provisional=True,
    )
    max_terminal_offset_bp: PolicyValue[int] = PolicyValue(
        value=25,
        justification="Chromosome ends terminate directly with canonical telomeric repeats; empirically adjusted as up to 16 bp of true telomere were not capped by exact hexamers (e.g. clipped leader sequence), so <= 25 bp tolerates this while excluding interstitial noise; BIO-001",
        provenance="@pgonzale60 biological calibration on reference benchmark (BIO-001)",
        provisional=True,
    )
    max_internal_gap_bp: PolicyValue[int] = PolicyValue(
        value=25,
        justification="Canonical telomeric repeat arrays are dense and contiguous; empirically adjusted as internal sequencing indels created gaps up to 12 bp within true telomeric arrays, so <= 25 bp preserves array continuity while breaking on non-telomeric sequence; BIO-001",
        provenance="@pgonzale60 biological calibration on reference benchmark (BIO-001)",
        provisional=True,
    )


@dataclass(frozen=True)
class StoichiometryPolicy:
    min_stem_core_ratio: PolicyValue[float] = PolicyValue(
        value=0.50,
        justification="Minimum stem read depth relative to core coverage to qualify as a substantial somatic stem; BIO-002",
        provenance="@pgonzale60 biological ratification on Issue #57; BIO-002",
        provisional=False,
    )
    min_germline_depth_ratio: PolicyValue[float] = PolicyValue(
        value=0.01,
        justification="Lower bound for germline-retained depth relative to somatic stem depth in mixed worm tissue; permits inclusion of lowest-coverage GRS (down to 1% of core coverage); BIO-002",
        provenance="@pgonzale60 biological ratification on Issue #57 (accepts down to 1% core coverage to include lowest-coverage GRS); BIO-002",
        provisional=False,
    )
    max_germline_depth_ratio: PolicyValue[float] = PolicyValue(
        value=0.33,
        justification="Upper bound for germline-retained depth relative to somatic stem depth in mixed worm tissue (0.33x core somatic coverage C_core); ratified on Issue #26.",
        provenance="Empirical reference benchmark data calibration ratified on Issue #26; BIO-002",
        provisional=False,
    )


@dataclass(frozen=True)
class PdePolicy:
    max_somatic_cap_hops: PolicyValue[int] = PolicyValue(
        value=2,
        justification="Maximum traversal hops from a somatic bifurcation port to identify an authentic terminal telomeric cap, accommodating short linker or bubble segments before the telomere tract; BIO-003",
        provenance="@pgonzale60 biological ratification on Issue #57; BIO-003",
        provisional=False,
    )


@dataclass(frozen=True)
class OntBridgingPolicy:
    min_anchor_bp: PolicyValue[int] = PolicyValue(
        value=1000,
        justification="Minimum alignment anchor length on flanking unitigs to ensure unambiguous physical anchoring across the satellite gap.",
        provenance="@pgonzale60 biological ratification on Issue #57; BIO-004",
        provisional=False,
    )
    max_boundary_slack_bp: PolicyValue[int] = PolicyValue(
        value=1500,
        justification="Maximum unaligned bases allowed between alignment terminus and the designated gap port boundary, accommodating repeat boundary fuzziness.",
        provenance="@pgonzale60 biological ratification on Issue #57; BIO-004",
        provisional=False,
    )
    min_mapq: PolicyValue[int] = PolicyValue(
        value=0,
        justification="Minimum mapping quality for repeat-flanking alignments; 0 permits repetitive satellite alignments when anchored by unique dovetail.",
        provenance="@pgonzale60 biological ratification on Issue #57",
        provisional=False,
    )
    min_spanning_reads: PolicyValue[int] = PolicyValue(
        value=2,
        justification="Minimum distinct spanning reads required to confirm physical bridging across the satellite array and declare END_COMPLETE.",
        provenance="@pgonzale60 biological ratification on Issue #57; BIO-004",
        provisional=False,
    )
    quantile_num: PolicyValue[int] = PolicyValue(
        value=4,
        justification="Number of quantiles for quartile and interquartile range (IQR) computation (quartiles = 4).",
        provenance="Mathematical definition of quartile distribution statistics",
        provisional=False,
    )


@dataclass(frozen=True)
class RdnaPolicy:
    min_45s_coverage_fraction: PolicyValue[float] = PolicyValue(
        value=0.15,
        justification="Minimum fraction of segment length covered by defined 45S rRNA regions (18S, 5.8S, 28S and intervening ITS) to classify as an authentic rDNA repeat cassette unitig; in benchmark assemblies, non-rDNA chromosome arms have <12.5% coverage (boundary hits <5%), whereas genuine rDNA unitigs and cassettes carry >=15% (typically 70-98%) 45S region sequence.",
        provenance="Biological calibration on reference benchmark Barrnap annotations (BIO-005)",
        provisional=False,
    )
    default_max_its_separation_bp: PolicyValue[int] = PolicyValue(
        value=5000,
        justification="Default fallback maximum unannotated separation distance between consecutive Barrnap 45S rRNA annotations (18S, 5.8S, 28S) to cluster into a continuous 45S rRNA region when empirical ITS size distribution cannot be inferred from the data (e.g. isolated hit).",
        provenance="Biological calibration on eukaryotic 45S rDNA ITS1/ITS2/IGS spacer architecture (BIO-005)",
        provisional=False,
    )
    min_its_threshold_floor_bp: PolicyValue[int] = PolicyValue(
        value=3000,
        justification="Minimum floor for empirical ITS clustering threshold distance to accommodate standard eukaryotic ITS1/ITS2 and short IGS spacers.",
        provenance="Biological calibration on eukaryotic 45S rDNA architecture (BIO-005)",
        provisional=False,
    )
    max_candidate_its_spacer_bp: PolicyValue[int] = PolicyValue(
        value=10000,
        justification="Upper bound for candidate within-array ITS/IGS spacers when distinguishing local operon spacing from chromosomal inter-locus genomic gaps.",
        provenance="Biological calibration on eukaryotic 45S rDNA architecture (BIO-005)",
        provisional=False,
    )
    its_separation_jump_ratio: PolicyValue[float] = PolicyValue(
        value=3.0,
        justification="Minimum ratio jump between sorted spacer sizes to detect the boundary between intra-array ITS spacers and inter-locus genomic gaps.",
        provenance="Empirical gap detection calibrated on benchmark assemblies (BIO-005)",
        provisional=False,
    )
    its_threshold_multiplier: PolicyValue[float] = PolicyValue(
        value=1.5,
        justification="Safety multiplier applied to maximum observed intra-array ITS spacer to define the clustering threshold.",
        provenance="Biological calibration on eukaryotic 45S rDNA architecture (BIO-005)",
        provisional=False,
    )


@dataclass(frozen=True)
class Policy:
    telomere: TelomerePolicy = TelomerePolicy()
    stoichiometry: StoichiometryPolicy = StoichiometryPolicy()
    pde: PdePolicy = PdePolicy()
    ont: OntBridgingPolicy = OntBridgingPolicy()
    rdna: RdnaPolicy = RdnaPolicy()


DEFAULT_POLICY = Policy()
