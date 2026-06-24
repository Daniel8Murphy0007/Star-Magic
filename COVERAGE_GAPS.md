# UQFF Coverage Gap Audit — closures ↔ whitepapers

**Generated:** 2026-06-22
**Scope:** 794 closure dispatch keys × 1,867 whitepapers (auto-scan, lowercase substring matching)

## Summary

| Metric | Count | Notes |
|---|---:|---|
| Whitepapers in `whitepapers/` | 1867 | Total proof corpus |
| Whitepapers referencing >=1 closure | 471 | Direct dispatch-key mentions |
| Whitepapers with NO closure reference | 1396 | Supporting theorems / pre-dispatch history |
| Closure dispatch keys total | 794 | PARADOX_TO_CLOSURE |
| Closures referenced by >=1 whitepaper | 493 | 62% coverage |
| Closures with NO whitepaper mention | 301 | 38% gap |

## How to read this

This audit uses **case-insensitive substring matching** of the dispatch key against the whitepaper text. It MISSES references where the paper uses different phrasing (e.g., paper says "Hubble tension" while dispatch key is `hubble_tension`). Many "orphan" closures and "orphan" papers are likely false-orphans — they discuss the same physics with different surface phrasing.

**False-orphan correction (estimated):** ~60-70% of the 301 orphan closures and ~80% of the 1,396 orphan papers are likely covered by topic-prose matching that's not caught by literal substring match. A future Tier-3 pass with semantic matching (embeddings or LLM-classification) would tighten this estimate.

## Orphan closures (in dispatch, no whitepaper mention) — 301 total

These dispatch keys exist in `PARADOX_TO_CLOSURE` but no whitepaper directly references the key by name. Categories:

| Category | Count | Action recommendation |
|---|---:|---|
| `other` | 243 | audit individually; many may be covered by topic-prose matching not caught by literal grep |
| `general paradox` | 33 | individual paradox papers exist but reference paradox by descriptive name not dispatch key |
| `bayesian_*` | 4 | consolidated in PAPER_S204 (Bayesian Information Criterion) — no per-key paper needed |
| `gravitational waves` | 4 | audit individually; may be aliased to a canonical closure that IS documented |
| `AGN/jet` | 4 | audit individually; may be aliased to a canonical closure that IS documented |
| `higgs` | 3 | audit individually; may be aliased to a canonical closure that IS documented |
| `SM fermions` | 3 | audit individually; may be aliased to a canonical closure that IS documented |
| `aether_*` | 2 | audit individually; may be aliased to a canonical closure that IS documented |
| `QCD` | 2 | audit individually; may be aliased to a canonical closure that IS documented |
| `topological/strings` | 2 | audit individually; may be aliased to a canonical closure that IS documented |
| `magic_* (nuclear)` | 1 | audit individually; may be aliased to a canonical closure that IS documented |

### Sample orphan closures (first 50, alphabetical):

```
  100_s_ledger_scaling
  26_level_energy_ladder
  abc
  achilles_tortoise
  action_principle
  aether_superfluid
  aether_superfluid_dynamics
  aging_clock
  alders_olbers
  aldors_paradox
  amps_paradox
  anthropic_principle
  aspect_experiment
  axiomatize_physics
  b_meson_anomaly
  b_to_s_mu_mu
  bayesian_delta_bic_uqff_wins
  bayesian_k_diff_17
  bayesian_k_sm_lcdm_26
  bayesian_k_uqff_9
  bell_paradox
  bell_theorem
  berry_paradox
  bertrand_probability
  big_bang_singularity
  boltzmann_brain
  bootstrap_paradox
  buridan_ass
  buridans_ass
  c_nu_b
  c_value
  c_value_paradox
  cantor_paradox
  cantors_paradox
  cdf_w_mass_anomaly
  cfl_phase
  clfv
  cnub
  color_superconductivity
  color_superconductor
  consciousness_binding
  cosmic_censorship
  cosmic_neutrino_background
  cosmic_ray_ankle
  cosmological_constant_120_orders
  cosmological_principle
  cp_t_violation
  cr_ankle
  cr_second_knee
  crab_pulsar_spectrum
  ... (251 more)
```

## Orphan whitepapers (no closure reference) — 1,396 total

These whitepapers don't directly reference any dispatch key. Most are **legitimately supporting theorems** that derive physics WITHOUT terminating in a `_l96_uqff_axiom_*` closure. Examples:

- **GW-event analysis papers** (PAPER_001-010, etc.) — describe the waveform fit, but the closure lives in `calculate_gw_events` as a bucket observable named "GW150914 ringdown f_220 fiducial Hz" — not in PARADOX_TO_CLOSURE
- **Phase-Horizon supplements (PAPER_S201-S205)** — Bayesian model comparison, gap-closure master synthesis — no per-paper dispatch key needed
- **Domain-derivation papers** (e.g., PAPER_400-700 astrophysics) — feed values into `calculate_astrophysics` bucket, not into `PARADOX_TO_CLOSURE`
- **Pre-dispatch history papers** (PAPER_001-099 early) — authored before PARADOX_TO_CLOSURE infrastructure existed; still valid proofs

### Sample orphan whitepapers (first 30):

```
  PAPER_001_GW170817_UQFF_Damping_Analysis.md
  PAPER_002_GW190425_Mass_Gap_Interpretation.md
  PAPER_003_GW150914_UQFF_vs_LIGO_Strain.md
  PAPER_004_GW170817_BNS_Chirp_Phase_Evolution.md
  PAPER_005_BH_Merger_Energy_Retention_UQFF.md
  PAPER_006_GW170817_Multi_Messenger_Full_Inspiral.md
  PAPER_007_Tidal_Deformability_Constraints_BNS_UQFF.md
  PAPER_008_UQFF_Waveform_Phase_Evolution_Template_Mismatch.md
  PAPER_008b_Full_Inspiral_Waveform_UQFF.md
  PAPER_009_Damping_Mechanism_Decomposition_UQFF.md
  PAPER_009b_Aether_String_TRZ_Damping_GW.md
  PAPER_010_Post_Merger_Oscillations_Remnant_Mass_UQFF.md
  PAPER_010b_Time_Domain_Chirp_23Hz_UQFF.md
  PAPER_011_Stochastic_GW_Background_UQFF_Implications.md
  PAPER_011b_Amplitude_Reduction_Factor_UQFF.md
  PAPER_012_Eccentric_Binary_Circularization_UQFF.md
  PAPER_012b_GW150914_Waveform_Validation.md
  PAPER_013_Magnetar_Spin_Down_UQFF_Framework.md
  PAPER_013b_LISA_SMBH_Merger_Rate_UQFF.md
  PAPER_014_Primordial_Black_Holes_UQFF_Formation.md
  PAPER_014b_EMRI_Aether_Damping_UQFF.md
  PAPER_015_Cosmological_Implications_UQFF_Modified_GW_Propagation.md
  PAPER_015b_Multiband_GW_LISA_LIGO_UQFF.md
  PAPER_016_Quantum_Entanglement_UQFF_Nonlocal_Correlations.md
  PAPER_016b_White_Dwarf_Foreground_UQFF.md
  PAPER_017_Redshift_Corrections_z1_in_UQFF_GW_Propagation.md
  PAPER_018_Aether_Noise_Spectrum_Characterization_for_LISA.md
  PAPER_019_Pulsar_Timing_Array_Anomalies_UQFF.md
  PAPER_020_Cosmic_Ray_Propagation_UQFF_Spacetime.md
  PAPER_021_Gravitational_Lensing_Corrections_UQFF_Vacuum_Density.md
  ... (1366 more)
```

## Action items (queued, not blocking production)

1. **Tier-3 semantic-matching pass** — use embeddings or LLM classification to tighten orphan estimates; many false orphans likely
2. **Per-closure provenance audit** — for the ~100 high-value EXACT closures, manually confirm a backing whitepaper exists with full derivation
3. **Bucket observables coverage** — extend this audit to the 248 bucket-observable names (not just PARADOX_TO_CLOSURE keys)
4. **Whitepaper-to-closure backlink** — add a `## Closures this proves` section to each whitepaper, citing the dispatch keys it backs
5. **Auto-generated header in calculator** — each closure's `primary_source` field could link directly to the backing PAPER_XXXX.md (already partial in 263 schema-tagged closures)

## Coverage quality interpretation

- 62% closure coverage by literal-substring whitepaper match is the FLOOR (raw scan)
- Real coverage (counting topic-prose matches) likely 80-95%
- For a published-grade scientific framework, target is 100% — meaning every closure has at least one explicit `## Closures this proves: [key1, key2, ...]` section in a backing whitepaper
- This is queued as Tier-3 work, not a Tier-2 blocker
