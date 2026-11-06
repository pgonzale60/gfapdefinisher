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
        provenance="@pgonzale60 biological calibration from somatic chromosome depth vs GRS; BIO-002",
        provisional=True,
    )
    min_germline_depth_ratio: PolicyValue[float] = PolicyValue(
        value=0.04,
        justification="Lower bound for germline-retained depth relative to somatic stem depth in mixed worm tissue; ~5% of cells are germline in unstaged culture, avoiding filtering out true GRS.",
        provenance="@pgonzale60 biological calibration from unstaged nematode cell proportions and GRS sensitivity (BIO-002)",
        provisional=True,
    )
    max_germline_depth_ratio: PolicyValue[float] = PolicyValue(
        value=0.40,
        justification="Upper bound for germline-retained depth relative to somatic stem depth in mixed worm tissue.",
        provenance="ha_analyze_count k-mer peak calibration on reference benchmark; BIO-002",
        provisional=True,
    )


@dataclass(frozen=True)
class PdePolicy:
    max_somatic_cap_hops: PolicyValue[int] = PolicyValue(
        value=2,
        justification="Maximum traversal hops from a somatic bifurcation port to identify an authentic terminal telomeric cap, accommodating short linker or bubble segments before the telomere tract; BIO-003",
        provenance="Empirical benchmark observation across curated nematode PDE junctions; BIO-003",
        provisional=True,
    )


@dataclass(frozen=True)
class Policy:
    telomere: TelomerePolicy = TelomerePolicy()
    stoichiometry: StoichiometryPolicy = StoichiometryPolicy()
    pde: PdePolicy = PdePolicy()


DEFAULT_POLICY = Policy()
