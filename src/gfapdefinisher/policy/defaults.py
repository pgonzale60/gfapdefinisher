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
    min_terminal_repeat_fraction: PolicyValue[float] = PolicyValue(
        value=0.60,
        justification="Minimum fraction of window composed of telomeric motifs to qualify as an authentic cap.",
        provenance="Empirical observation across curated nematode T2T assemblies",
        provisional=True,
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
        value=10,
        justification="Chromosome ends terminate directly with canonical telomeric repeats; non-repeat leader sequence must be <= 10 bp to distinguish genuine terminal caps from subtelomeric/interstitial repeats and tolerate only minimal sequencing or clipping noise; BIO-001",
        provenance="@pgonzale60 biological calibration on reference benchmark (BIO-001)",
        provisional=True,
    )
    max_internal_gap_bp: PolicyValue[int] = PolicyValue(
        value=10,
        justification="Canonical telomeric repeat arrays are dense and contiguous; non-repeat gaps between adjacent repeat units exceeding 10 bp indicate interrupting non-telomeric sequence rather than minor sequencing error; BIO-001",
        provenance="@pgonzale60 biological calibration on reference benchmark (BIO-001)",
        provisional=True,
    )


@dataclass(frozen=True)
class StoichiometryPolicy:
    min_germline_depth_ratio: PolicyValue[float] = PolicyValue(
        value=0.15,
        justification="Lower bound for germline-retained depth relative to somatic stem depth in mixed worm tissue.",
        provenance="ha_analyze_count k-mer peak calibration on reference benchmark; BIO-002",
        provisional=True,
    )
    max_germline_depth_ratio: PolicyValue[float] = PolicyValue(
        value=0.45,
        justification="Upper bound for germline-retained depth relative to somatic stem depth in mixed worm tissue.",
        provenance="ha_analyze_count k-mer peak calibration on reference benchmark; BIO-002",
        provisional=True,
    )


@dataclass(frozen=True)
class Policy:
    telomere: TelomerePolicy = TelomerePolicy()
    stoichiometry: StoichiometryPolicy = StoichiometryPolicy()


DEFAULT_POLICY = Policy()
