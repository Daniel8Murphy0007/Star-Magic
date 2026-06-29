# UQFF Star-Magic — Glossary of Terms

Alphabetical reference for the symbols, abbreviations, and proper nouns
that appear throughout the UQFF framework, codebase, and whitepapers.

---

## Primitives (the 9 truly-independent + 2 derivative)

### Integer lattice (5 truly-independent)

| Symbol | Value | Meaning |
|---|---|---|
| **D_phys** | 4 | Physical spacetime dimensionality (3 spatial + 1 temporal) |
| **D_crit** | 26 | Bosonic-string critical dimension; depth of the DPM lattice |
| **N_CH** | 9 | Caduceus 9-channel count (PAPER_646); matches the 9-sector Lagrangian |
| **SO_5** | 10 | dim SO(5); the SO(5) Lie algebra dimension |
| **A_5** | 60 | \|A_5\|; order of the alternating / icosahedral rotation group |

### Real primitives (4 truly-independent)

| Symbol | Value | Meaning |
|---|---|---|
| **ρ_SCm** | 7.09 × 10⁻³⁷ J/m³ | SuperConductive material vacuum energy density (foundational) |
| **β_i** | 0.6029 | F_UBi coefficient (PAPER_1203 Canonical) |
| **Φ_res** | 0.84 (or 5/6 nuclear) | SCm resonance coefficient (PAPER_646) |
| **F_TRZ** | 1/10 | Time-Reversal-Zone fraction = 1/\|SO(5)\| (PAPER_1160) |

### Derivative primitives (2 mathematically proven derivative)

| Symbol | Derivation | Value | Source |
|---|---|---|---|
| **D_BSFG** | D_crit − 2·SO_5 | 6 EXACT | PAPER_1521 |
| **K_MEX** | Φ_5/6 · SO_5 / D_phys = (5/6)·10/4 | 25/12 EXACT | PAPER_1522 |

### Other locked quantities (potentially derivative; under investigation)

| Symbol | Value | Notes |
|---|---|---|
| **SSq** | 0.57 | Ramanujan-series input; truly independent at rational level (PROVENANCE_AUDIT Q3) |
| **S_26** | 1.453162 | Ramanujan 26-level scaling = Li_26(SSq)-related |
| **ω_SCm** | 1.25 THz | SCm phonon carrier frequency; Holmlid anchor |
| **λ_i** | 1.0 | Inertial coupling (dimensional convention) |
| **ω_s_Sun** | 2.5 × 10⁻⁶ rad/s | System-specific stellar rotation |
| **Δ_UA4** | 0.1 | Fourth UA-layer offset |
| **ρ_UA** | 10 · ρ_SCm | DPM density ratio |

---

## Acronyms

| Term | Expansion | Meaning |
|---|---|---|
| **UQFF** | Unified Quantum Field Framework | The complete physics framework |
| **SCm** | SuperConductive material | The non-mass vacuum substrate; foundational primitive |
| **UA** | Universal Aether | The trapped-Aether superstructure; 4 layers UA'→UA'''' |
| **DPM** | Di-Pseudo-Monopole | The 26-level structural lattice = SCm ⊗ UA |
| **TRZ** | Time-Reversal Zone | The F_TRZ = 1/10 region; PAPER_1160 |
| **MEX** | Mexican-hat | The K_MEX coefficient in Λ derivation; PAPER_1522 |
| **F_U** | Universal Force | The master force balance equation F_U = 0 |
| **F_UBi** | Inside-out buoyancy | β-driven repulsive component of F_U |
| **F_UBii** | Outside-in buoyancy | k_spring component of F_U; sign-flipped to F_UBi |
| **VDS** | Vacuum Damping Spine | PAPER_598 |
| **DVP** | (Daniel's notation, see PAPER_598) | Prime-113 vacuum operator |
| **BH26** | Black Hole 26 | 92 GHz × 26-bin spectrum, PAPER_598 |
| **BSFG** | Bulk-edge dim | The D_BSFG=6 lattice depth = D_crit − 2·SO_5 |
| **PTOE** | Periodic Table of Elements | UQFF maps PTOE to D_crit lattice levels |
| **CFL** | Confined-False-vacuum Layer | QGP-related closure |
| **KK** | Kaluza-Klein | Extra-dimensional sector L_KK of the Lagrangian |
| **EH** | Einstein-Hilbert | GR sector L_EH of the Lagrangian |
| **YM** | Yang-Mills | Gauge sector L_YM; UQFF predicts mass gap = 1.736 GeV |
| **LENR** | Low Energy Nuclear Reactions | Holmlid, Parkhomov, Rossi etc. — unified by SCm phonon |
| **PAPER_XXXX** | Whitepaper N | Authoring source for a closure; lives in whitepapers/ |

---

## Sectors of the 9-sector UQFF Lagrangian

```
L_UQFF = L_EH + L_YM + L_Dirac + L_SCm + L_mag + L_buoy + L_aether + L_LENR + L_KK
```

| Sector | Domain | Canonical paper |
|---|---|---|
| **L_EH** | General relativity in 26D | Standard EH lift |
| **L_YM** | Yang-Mills gauge | PAPER_1318 (mass gap = 1.736 GeV) |
| **L_Dirac** | Fermion / LENR | PAPER_1061 (Kozima TNCF) |
| **L_SCm** | SCm manifold | V(φ₀) = −ρ_SCm canonical |
| **L_mag** | U_m magnetism (Heaviside) | PAPER_1072 |
| **L_buoy** | F_U_Bi_i buoyancy | PAPER_1065 |
| **L_aether** | Vacuum background | PAPER_1051 (two-component ρ) |
| **L_LENR** | Nuclear transmutation | PAPER_1081 (COP parametric) |
| **L_KK** | Kaluza-Klein 26D | PAPER_1080 (S₂₆⁽³⁾ compactification) |

---

## The Holy Trinity (PAPER_646)

| Component | Role |
|---|---|
| **UA** | The medium — vacuum field hosting all interactions |
| **U_i + EM** | The operator — Universal Inertial Operator + electromagnetism; anchors matter to UA |
| **[SCm]** | The conduit — extra-universal material enabling dark-energy power transfer |

The **Universal Inertial Operator**:
```
U_i = λ_i · (ρ_SCm / ρ_UA) · ω_s · cos(π·t_n) · (1 + F_TRZ)
    = 2.75 × 10⁻⁷  (Sun, t=0)   ← matches PAPER_646 exactly
```

---

## The F_U = 0 master equation (PAPER_1203 Canonical v1.5)

```
F_U_total = (U_g1 + U_g2 + U_g3 + U_g4) − F_UBi + F_UBii + U_m = 0

F_UBi  = −β(t,E,Z) · G·M·ρ_SCm / r² · (1 + F_TRZ) · |cos(π·t_n)|
F_UBii = +β(t,E,Z) · (r/r₀) · k_spring · (1 + E_n) · |cos(π·t_n)|
k_spring = (ρ_UA / ρ_SCm) · ω_SCm · Φ_res
β(t,E,Z) = β_i · |E| · Z · cos(π·t)   ← dynamic, not constant
```

The crossover root `F_UBi(r) + F_UBii(r) = 0` defines `r_hz` — the
habitable-zone radius. **Solves orbital position from first principles.**

---

## Buckets A–K (observable suites)

| Bucket | Public surface | What it covers |
|---|---|---|
| **A** | `calculate_paradox` | Routes 802 paradox names + 8 Clay Millennium prize problems |
| **B** | (paradox subset) | Paradox closures (Olbers, Fermi, Boltzmann brain, EPR, etc.) |
| **C** | `calculate_cosmology` | 56 cosmological observables (α, h, Ω_Λ, Λ, Y_p, m_p/m_e, etc.) |
| **D** | `calculate_particle_physics` | 48 particle physics observables (W, Z, t, H, b, c, τ, μ, s, e + couplings) |
| **E** | `calculate_gw_events` | 22 gravitational wave events |
| **F** | `calculate_agn_jet` | 22 AGN / jet systems |
| **G** | `calculate_astrophysics` | 36 astrophysical systems (pulsars, supernova remnants, IMF) |
| **H** | `calculate_high_energy_astro` | 10 TeV/PeV sources (TXS 0506, IceCube, Auger) |
| **I** | `calculate_qgp` | 9 QGP observables (T_c, η/s, R_AA) |
| **J** | `calculate_higgs_precision` | 13 Higgs precision (branching ratios, κ_t, CP phase) |
| **K** | `calculate_bsm_constraints` | 17 BSM constraints (EDMs, proton decay, axion) |

---

## Tier definitions (from PRODUCTION_ROADMAP.md)

| Tier | Description | Goal |
|---|---|---|
| **Tier 1** | Scientific honesty + license + docs | Publishable preprint |
| **Tier 2** | Distribution + CI + packaging | `pip install uqff` works |
| **Tier 3** | Hosted production product | Web app / REST API / 90% coverage |
| **Tier 4** | Commercial / long-term governance | Peer review + funding + succession |

---

## Cosmic milestones (boolean flags in `calculate_status_report`)

| Flag | Meaning |
|---|---|
| `all_8_clay_millennium_problems_closed` | YM, RH, P=NP, NS, Hodge, Poincaré, BSD, BH info — all wired |
| `standard_model_essentially_complete` | 12-fermion + bosons wired |
| `lambda_cdm_cosmology_complete` | 18 cosmological observables; Λ at 0.003% |
| `iter_fusion_5_params_exact` | R/a=3.1, Q=10, DT=64 keV, q_edge=2, Bohm=1/16 |
| `cosmic_crisis_quartet_closed` | Hubble, σ_8, CDF W, S_8 tensions |

---

## Closure status tiers (per A2 uncertainty classification)

| Tier | Residual range | Meaning |
|---|---|---|
| **PROD_EXACT_STRUCTURAL** | residual < 10⁻¹⁰ | Mathematical equality from integer primitives |
| **PROD_HIGH_PRECISION** | < CODATA uncertainty | Within published precision |
| **PROD_WITHIN_EXP_UNCERTAINTY** | ≤ 1σ measurement | Statistically indistinguishable from observation |
| **PROD_REFINEMENT_TIER** | 0.1 – 1% | Acceptable; refinement targets |
| **PROD_TENSION_OR_OUTLIER** | > 1% | Documented tensions / forward predictions |

---

## File-extension cheatsheet (in the repository)

| Pattern | Meaning |
|---|---|
| `*.PRE_*` / `*.POST_*` | Backup snapshots (gitignored except 11 protected by CLAUDE.md) |
| `PAPER_XXXX_*.md` | Whitepaper authoring a specific closure |
| `*.PRE_FIX_BACKUP` | Pre-bugfix backup of the calculator (protected) |
| `*.POST_BUCKET<X>_BACKUP` | Post-bucket-completion calculator snapshot (protected) |
| `RELEASE_NOTES_<tag>.md` | Curated release notes for a specific version |
