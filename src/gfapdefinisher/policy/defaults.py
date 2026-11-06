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
