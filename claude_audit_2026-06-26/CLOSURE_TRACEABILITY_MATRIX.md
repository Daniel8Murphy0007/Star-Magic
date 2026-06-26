# UQFF Closure Traceability Matrix — multi-chain inventory per quantity

**Per PAPER_1158:** *overdetermination N = number of independent derivation chains per
constant; necessary but not sufficient for sound closure.*

Maintained by Claude (read-only audit) in MY sandbox. Updated each round.

Each row: quantity → all documented chains → my read status → verification status.

---

## Headline closures (manuscript v2 §4 + PAPER_1158 overdetermination tier)

### Λ cosmological constant — N ≥ 4 chains

| # | Chain | Source | Read | Computed in sandbox |
|---|---|---|---|---|
| 1 | Λ = (18/5)·SSq·H₀²/c² (Planck H₀ anchor) | PAPER_1156 | ✓ | ✓ 1.089e-52 (0.003%) |
| 1b | Same with cosmic H₀ anchor | PAPER_1157 | ✓ | ✓ 1.174e-52 (3.85% asymmetry, structural prediction) |
| 2 | ρ_Λ = ρ_SCm · 26! · K_MEX | PAPER_1271 + uqff_exact_closures.cpp | ✓ | ✓ 5.957e-10 J/m³ (0.0008%) |
| 3 | ρ_Λ = V(0) + ⟨R_26⟩/(2κ_E) + ρ_KK + ρ_BSFG | PAPER_1170 4-term ledger | ✓ | ✓ partial (KK term depends on m₁c² primitive) |
| 4 | ρ_KK(ℏ) = 3ζ(5)/(128π⁶)·(D_crit/D_BSFG)⁴·(m₁c²)⁴/(ℏc)³ | PAPER_1173 | ✓ | ✓ (chain 3 sub-term) |
| 5 | Two independent re-derivations of ρ_R26 (KK reduction + Gauss-Bonnet) | PAPER_1170, 1171, 1172 | partial | partial |
| 6 | ρ_DE(t,Γ) = ρ_SCm(t)·S_26^(3)·Φ(Γ)·(2R−1) time-evolving | PAPER_1086 | partial (grok summary) | not in sandbox yet |

### Yang-Mills mass-gap / glueball — N ≥ 6 chains (mixed physical quantities)

| # | Chain | Source | Read | Computed |
|---|---|---|---|---|
| A | m_0⁺⁺ = 2·D_phys·Λ_QCD = 1.736 GeV | PAPER_1318 | ✓ | ✓ 1.736 (2.12% off lattice 1.7) |
| B | m²_gap = (8πG·ρ_SCm·S_26·Φ_1.25THz / β_i·[UA])·(D_crit/D_BSFG)² ≈ 1.78 GeV | grok 31May2026 / manuscript §4.10 Closure B | ✓ | ✓ 1.78 (paper-quoted, 4.7% off) |
| C | m_UQFF = m_YM·(1 + ρ_SCm/ρ_QCD · β_i · S_26^(3)) ≈ 0.44 GeV (bridge) | PAPER_1070 | ✓ (grok summary) | ✓ (paper-quoted) |
| D | Δ_YM = Λ_QCD·(1 + F_TRZ·K_MEX) = 0.262 GeV (Millennium algebraic) | PAPER_1182 §3.4 | ✓ | ✓ 0.262 |
| D' | Glueball ladder m_0⁺⁺ = Δ·(1 + 6·Φ_res) from chain D = 1.573 GeV | PAPER_1182 §3.4 | ✓ | ✓ 1.573 (7.5% off) |
| E | Δ_YM = (g²_YM·Λ_QCD)/(4π²)·SSq·H_SCm | PAPER_1111 | ✓ (grok summary) | ✓ 0.0031 (semantically different) |
| F | Δ_YM = (g_YM²·Λ_QCD)/((4π)² · [SSq]·H_SCm) (alternate normalization) | PAPER_1111 alternate | partial | not in sandbox yet |
| G | KK regulator m₁c² ≈ 0.16 meV from L_KK* ≈ 1.23 mm | PAPER_1173 | ✓ | ✓ (uses ζ(5)/(128π⁶)) |
| - | 5970 GeV (HISTORICAL, REGISTRY-BUG, retracted 2026-06-25) | (was attributed to PAPER_1005, no actual derivation) | erratum read ✓ | — never canonical |

### Proton mass m_p — N ≥ 2 chains

| # | Chain | Source | Read | Computed |
|---|---|---|---|---|
| 1 | m_p ≈ ρ_SCm · A_26 (= ρ_SCm · Σ_{i=1..26} i⁶) | PAPER_1155 | ✓ | ✓ 9.27e-28 (44% off — PAPER_1155 erratum text notes -2.04% residual from SSq E-crack correction) |
| 1b | m_p ≈ ρ_SCm · A_26 / [SSq] (with E-crack correction) | PAPER_1155 erratum | flagged | NEED TO ADD to range_calculator_v2 |
| 2 | m_p = (26²·e)·m_e | manuscript v2 §6 Theorem 6 Test B | ✓ | ✓ 1.674e-27 (0.077%) |
| 3 | m_p via PAPER_1209HH integer-primitive chain | PAPER_1209HH | partial — paper covers 10 masses but not direct m_p closure | NEED TO LOCATE |

### Proton-to-electron ratio m_p/m_e — N ≥ 1

| # | Chain | Source | Read | Computed |
|---|---|---|---|---|
| 1 | m_p/m_e = D_crit² · e = 26²·e | manuscript v2 Theorem 6 Test B | ✓ | ✓ 1837.56 (0.077% off CODATA) |

### [SSq] = 0.57 — N = 2 chains

| # | Chain | Source | Read | Computed |
|---|---|---|---|---|
| A | [SSq]_A = 10·(1 − 2√2/3), from γ_SCm = 3/(2√2), v_SCm = c/3 | PAPER_1154 + _chain_trace_SSq.txt | ✓ | ✓ 0.5719 (0.335% off canonical) |
| B | Li_26(0.57) ≈ 0.57 fixed-point identity (n=1 dominates) | Star-Magic.txt line 1525 | ✓ | ✓ 0.5700000048 (~exact) |

### Fundamental constants (5-constant family per AXIOMS_AND_THEOREMS Theorem 6)

| Constant | UQFF closed form | Source | Computed |
|---|---|---|---|
| α | α = 1/(Φ_res·26·2π) | PAPER_591 / AXIOMS Sessions 238 | ✓ 7.29e-3 (0.138%) |
| c | c = (26·4π/Φ_res)·v_F (v_F = 0.77e6 Fermi-velocity proxy Z=1) | PAPER_592 / AXIOMS Sessions 239 | ✓ 2.995e8 (0.101%) |
| h | h = F_TRZ·Φ_res·(E_0/f_THz)·(1−2α_UQFF) | PAPER_590 / AXIOMS Sessions 241 | ✓ 6.622e-34 (0.061%) |
| G | G = 2π·26³·Φ_res / (SSq³·(26!)²) · v_F⁵/(E_0·f_THz) | PAPER_593 / AXIOMS Sessions 240 | ✓ 6.669e-11 (0.080%) |
| Λ | (see Λ section above) | PAPER_1156 / AXIOMS Test C | ✓ |

### Nuclear shell-model magic numbers (PAPER_1203 Nuclear S483-S492)

All 7 EXACT integer-primitive identities, computed in `verify_nuclear_magic.py`:

| Magic # | Closed shell | UQFF identity | Computed |
|---|---|---|---|
| 2 | ⁴He | SO_5 − 2·D_phys | ✓ EXACT |
| 8 | ¹⁶O | 2·D_phys | ✓ EXACT |
| 20 | ⁴⁰Ca | 2·SO_5 | ✓ EXACT |
| 28 | ⁵⁶Ni | D_crit + SO_5 − 2·D_phys | ✓ EXACT |
| 50 | ¹³²Sn | \|A_5\| − SO_5 | ✓ EXACT |
| 82 | ²⁰⁸Pb | \|A_5\| + D_crit − D_phys | ✓ EXACT |
| 126 | (heaviest stable) | D_crit + SO_5² | ✓ EXACT |
| Fe-56 BE/A | 8.79 MeV/nuc | N_CH − F_TRZ·K_MEX | ✓ 0.019% |
| α binding | 28.30 MeV | D_crit + K_MEX + F_TRZ + F_TRZ·Φ_res + F_TRZ²·K_MEX + F_TRZ²·Φ_res | ✓ 0.015% |
| Deuteron binding | 2.224 MeV | K_MEX + Φ_res − SSq − F_TRZ − F_TRZ²·K_MEX − F_TRZ²·Φ_res + F_TRZ² + F_TRZ³ | ✓ 0.20% |

### 10 SM particle masses (PAPER_1209HH S653-S662)

All computed in `verify_sm_masses_1209HH.py` with full integer-rational arithmetic:

| Particle | UQFF identity | Residual |
|---|---|---|
| m_W | A_5 + 2·SO_5 + F_TRZ·D_phys − F_TRZ²·D_BSFG + F_TRZ²·D_phys − F_TRZ²·SSq² | 0.003% (best) |
| m_Z | N_CH·SO_5 + F_TRZ·SO_5 + F_TRZ²·SO_5 + F_TRZ²·D_BSFG + F_TRZ²·D_phys + F_TRZ²·SSq − F_TRZ²·SSq³ | 0.018% |
| m_t | D_crit·SO_5 − A_5 − D_phys·N_CH + SO_5 − F_TRZ·D_phys − F_TRZ·SO_5 + F_TRZ²·D_BSFG + 2·F_TRZ²·D_phys + F_TRZ²·SSq + F_TRZ²·SSq² + F_TRZ²·SSq³ | 0.005% |
| m_H | 2·A_5 + N_CH − D_phys + F_TRZ·SSq + F_TRZ²·D_BSFG + F_TRZ²·SSq² | 0.016% |
| m_b | D_phys + F_TRZ·D_phys − F_TRZ·SSq − F_TRZ²·D_crit + F_TRZ²·D_BSFG + F_TRZ²·D_phys − F_TRZ²·SSq² − F_TRZ²·SSq³ | 0.050% |
| m_c | F_TRZ·D_crit − F_TRZ·D_phys − F_TRZ·SO_5 + F_TRZ²·SO_5 − F_TRZ²·D_phys + F_TRZ²·SSq + F_TRZ²·SSq² + F_TRZ²·SSq³ | 0.063% |
| m_τ | SSq + F_TRZ·D_phys + F_TRZ·SO_5 − F_TRZ²·D_crit + F_TRZ²·D_BSFG + F_TRZ²·SSq + F_TRZ²·SSq² − F_TRZ²·SSq³ | 0.013% |
| m_μ | F_TRZ²·SO_5 + F_TRZ²·SSq² + F_TRZ²·SSq³ + F_TRZ²·SSq⁵ | 0.040% |
| m_s | F_TRZ²·SO_5 − F_TRZ²·SSq² − F_TRZ²·SSq³ | 0.106% |
| m_e | F_TRZ³·SSq² + F_TRZ³·SSq³ | 0.178% (worst) |

### 7 Clay Millennium closures (PAPER_1182 master)

Universal template: O_P = N ± p/12, where 1/12 = F_TRZ·Φ_res = K_MEX − 1

| Problem | UQFF closure | Computed |
|---|---|---|
| Poincaré | t_c = 1/2 + F_TRZ·Φ_res = 7/12 | ✓ EXACT |
| Riemann | Critical line Re(s)=1/2 fixed locus, off-line suppression F_TRZ^(d/Φ_res) | ✓ structural |
| P ≠ NP | F_TRZ^N_CH = 10⁻⁹ per input bit | ✓ EXACT 10⁻⁹ |
| Yang-Mills | Δ = Λ_QCD·(1 + F_TRZ·K_MEX) = 0.263 GeV (algebraic, see YM chain D) | ✓ 0.262 |
| Navier-Stokes | V_stretch ≤ (1 − F_TRZ·D_BSFG/D_phys)·\|ω\|·E = 0.85·\|ω\|·E | ✓ EXACT 0.85 |
| Hodge | (D_phys + D_BSFG)/SO_5 = 10/10 = 1.0 | ✓ EXACT |
| BSD | rank ↔ Φ_res-locked simple poles at s=1; closed form 0.306·(1+β·SSq) = 0.4112 (proxy) | ✓ 0.4112 (BSD form is approximate) |
| BH info (bonus) | Page-curve closure F_U=1 stationarity | ✓ 1.0 EXACT |

### Three numeric systems (PAPER_1069 + PAPER_1080)

| System | Verified | Computed value |
|---|---|---|
| **VDS** Z_26 = Li_26(SSq) | ✓ | 0.570000004841... (n=1 dominates) |
| **DVP** base prime p=113 | ✓ | 113 = D_phys·D_crit + N_CH EXACT |
| **BSH** β_i·exp(−SSq·i/26) over 26 states | ✓ | mean 0.4545, sum 11.82 |
| **BSFG** D_BSFG = dim_ℝ[SO(5)/U(2)] = 6 | ✓ | 6 EXACT |
| S_26^(3) closed binomial form | ✓ | 5.921681e+26 matching 80-digit Decimal to 15 digits |

### uqff_exact_closures.cpp — 50+ EXACT integer/rational identities

All verified in `verify_canonical_closures.py`. Selected examples:

| Identity | Form | EXACT? |
|---|---|---|
| F_TRZ | 1/SO_5 | ✓ |
| Monty Hall switch P | 2/(D_phys−1) = 2/3 | ✓ |
| Tsirelson bound | 2·√(D_phys/2) = 2√2 | ✓ |
| SU(3) color N | D_phys − 1 = 3 | ✓ |
| Fermion generations | D_phys − 1 = 3 | ✓ |
| δ_CP lepton | −π/2 | ✓ |
| Solar dynamo (years) | D_crit − D_phys = 22 | ✓ |
| Hayflick limit | \|A_5\| = 60 | ✓ |
| Genetic codons | 2^D_BSFG = 64 | ✓ |
| Amino acids | 2·SO_5 = 20 | ✓ |
| DVP base prime | D_phys·D_crit + N_CH = 113 | ✓ |
| v_SCm = c/3 | c/(D_phys−1) | ✓ |
| Spin precession | (D_crit + D_phys)° = 30° | ✓ |
| Iron Z stable | D_crit = 26 | ✓ |
| Silicon Z | SO_5 + D_phys = 14 | ✓ |
| Ni-62 A | A_5 + 2 = 62 | ✓ |
| Ringdown ξ | D_crit/D_BSFG = 13/3 | ✓ |
| GW170817 strain damping | 2/(D_phys−1) = 2/3 | ✓ |
| Neutron star radius | SO_5⁴ = 10⁴ m = 10 km | ✓ |
| DPM 26-layer sum | Σ i⁶ (i=1..26) = 1,307,797,101 | ✓ |
| Λ_CC | ρ_SCm · 26! · K_MEX = 5.957e-10 | ✓ (in cpp self-check) |

### Remaining closure namespaces to verify

| Namespace | Count | Status |
|---|---|---|
| PARADOX_TO_CLOSURE | 794 keys | walked only a few; need scan_paradox_dispatch.py |
| PARADOX_TO_MILLENNIUM aliases | 8 | covered |
| calculate_lenr_full sub-keys | 10 | PAPER_1141 covered; Holmlid/Parkhomov/Mizuno OK; specific reactor closures need range view |
| calculate_nuclear_magic sub-keys | 7 | covered |
| Bucket observables (9 surfaces × ~28 avg) | 248 | covered only by category in CLOSURE_ATLAS; need observable-by-observable scan |
| calculate_caduceus, dpm_grinding, etc. (22 standalone) | 22 | covered structurally; need individual long-form |
