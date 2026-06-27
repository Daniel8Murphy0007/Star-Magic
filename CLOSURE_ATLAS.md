# UQFF Closure Atlas — Where everything lives

**Date generated:** 2026-06-22 (v5.27.2)
**Author:** Daniel T. Murphy + closure-mapping pass
**Purpose:** Single authoritative reference for the location of every closure,
proof, derivation, and theorem in the UQFF corpus.

---

## TL;DR

```
   794  paradox dispatch keys (PARADOX_TO_CLOSURE)
+    8  Clay Millennium derivation functions (_MILLENNIUM_DERIVE)
+  248  bucket observables (across 9 surfaces)
+   22  standalone deep-content surfaces
+1,867  whitepapers (axiom proof sets + supporting theorems)
+  368  C++ reference closures
+  857  fidelity gate tests
─────
≈4,164 distinct named closure/proof artifacts
```

---

## 1. The 8 Clay Millennium Prize Problems

**File:** `uqff_pure_calculator.py`
**Dicts:** `PARADOX_TO_MILLENNIUM` (alias map) + `_MILLENNIUM_DERIVE` (functions)

| Public paradox key | Routes to | Derivation function | UQFF value |
|---|---|---|---|
| `yang_mills_mass_gap` | yang_mills | `_millennium_yang_mills_derive` | **1.736 GeV** |
| `riemann_hypothesis` | riemann | `_millennium_riemann_derive` | t_10000 EXACT (9877.78265) |
| `bsd_conjecture` | bsd | `_millennium_bsd_derive` | 0.30598 (Cremona 37a1, 0.005%) |
| `navier_stokes_smoothness` | navier_stokes | `_millennium_navier_stokes_derive` | enstrophy cap 0.85 |
| `hodge_conjecture` | hodge | `_millennium_hodge_derive` | 1.0 EXACT |
| `poincare_conjecture` | poincare | `_millennium_poincare_derive` | 7/12 |
| `p_vs_np` | p_vs_np | `_millennium_p_vs_np_derive` | 1 − 10⁻⁹ |
| `info_paradox` | black_hole_info | `_millennium_black_hole_info_derive` | 0.99596 |

**Callable via:**

```bash
uqff predict yang_mills
uqff predict riemann
uqff predict bsd
... etc
```

Or in Python:

```python
u.calculate_paradox({"paradox": "yang_mills"})
```

**Backing whitepapers:** PAPER_1005 (YM), PAPER_1182 (Riemann), PAPER_1183 (paradoxes),
PAPER_1230 (Hodge), and ~30 more across `whitepapers/PAPER_1200_*` and `PAPER_1300_*`.

---

## 2. Paradox dispatch table (PARADOX_TO_CLOSURE)

**File:** `uqff_pure_calculator.py` lines ~39,000-48,000

- 794 dispatch keys → 616 unique closure functions (178 aliases)
- Each closure is `_l96_uqff_axiom_<name>_closure() -> dict`
- All keys are **lowercase** (case-sensitivity bug otherwise)
- The dispatcher (`_paradox_proof`) normalizes via `.lower().strip()`

**Browse:**

```bash
uqff list                       # the 794 paradox dispatch keys
uqff list --filter <substr>     # filter
uqff list --all                 # include all 5 namespaces (~1,080 names)
```

**Discoverable categories within PARADOX_TO_CLOSURE:**

| Category | Approx count | Examples |
|---|---:|---|
| Quantum mechanics paradoxes | ~80 | EPR, GHZ, Hardy, Kochen-Specker, Bell, Wigner-friend |
| Relativity / black-hole paradoxes | ~50 | firewall, amps, page_curve, hawking_info |
| Cosmology paradoxes | ~40 | hubble_tension, sigma_8, cdf_w, dark_flow, Λ_problem |
| Quantum field theory | ~60 | strong_cp_problem, axion, gauge_coupling, vacuum_catastrophe |
| Standard Model identities | ~80 | m_W, m_Z, m_t, fermion masses, neutrino splittings |
| Nuclear / atomic | ~50 | magic numbers, BE/A, alpha, rydberg, proton_radius |
| LENR / fusion | ~30 | star_magic_reactor_resonance, ker_chain, transmut |
| Astrophysics | ~80 | psr_j0030, crab, eta_car, supernova, IMF |
| Statistical / information | ~30 | monty_hall, sleeping_beauty, simpson, bertrand |
| Mathematical | ~50 | hilbert_hotel, banach_tarski, peano, inverse_galois |
| Engineering / fusion / quantum-computing | ~40 | iter_*, surface_code, qubit_threshold |
| Other paradoxes | ~204 | misc (incl. PAPER_S### supplements) |

---

## 3. The 9 bucket observable surfaces (248 named observables)

**File:** `uqff_pure_calculator.py` (each `calculate_X` returns `{'value': {'observables': [...]}}`)

| Surface | Observables | Domain |
|---|---:|---|
| `calculate_cosmology` | 60 | α, h, Ω_Λ, Λ, m_p/m_e, Y_p, τ_reion, σ_8, ... |
| `calculate_particle_physics` | 49 | m_W, m_Z, m_top, m_H, CKM, neutrino splittings, ... |
| `calculate_astrophysics` | 43 | PSR J0030 mass, Crab spin-down, Eta Carinae, ... |
| `calculate_agn_jet` | 23 | 3C273 Eddington, TON618 mass, SgrA* shadow, M-σ, ... |
| `calculate_gw_events` | 22 | GW150914, GW170817, NANOGrav 15yr, LIGO O5, ... |
| `calculate_bsm_constraints` | 19 | electron/neutron EDM, proton decay, axion, ... |
| `calculate_higgs_precision` | 13 | H→γγ, H→ZZ, H→WW branching, κ_t, CP phase |
| `calculate_high_energy_astro` | 10 | TXS 0506+056, IceCube ν, Auger dipole, Amaterasu |
| `calculate_qgp` | 9 | T_c, η/s, ALICE PbPb dN_ch/dη |

**Discoverable via:**

```bash
uqff search <observable_substring>           # finds across all 9 surfaces
python -c "u.calculate_cosmology({})['value']['observables']"   # full list
```

Each observable record: `{observable, uqff_derived, anchor, residual_pct}`.

---

## 4. 22 standalone deep-content surfaces

Each returns structured dicts with named sub-fields, not paradox dispatch.

| Surface | Sub-keys |
|---|---|
| `calculate_lenr_full` | holmlid_D_minus_1, parkhomov_NiH, pons_fleischmann_PdD, mizuno_NiD, rossi_early/x/sk, **star_magic_reactor**, widom_larsen, kozima_tncf |
| `calculate_lenr` | ker_chain, parkhomov_1hr_W, pons_fleischmann_W |
| `calculate_nuclear_magic` | magic_numbers (7), magic_expected, magic_exact, magic_all_exact, be_per_a_peak_MeV, alpha_binding_MeV, deuteron, ... |
| `calculate_caduceus` | 26 pinch_points, pi_decimal_sequence (8 digits), sum_amplitudes |
| `calculate_dpm_grinding` | 5 steps of SCm→UA→UA''''' grinding |
| `calculate_ua_layers` | UA_prime, UA_double, UA_triple, UA_quad, F_DPM |
| `calculate_universal_inertial_operator` | U_i_dimensionless (= 2.75e-7 Sun), U_i_product_J_per_m3 |
| `calculate_vacuum_ledger` | V0, R26_term, ρ_KK, ρ_BSFG → Λ |
| `calculate_f_u_zero` | F_UBi, F_UBii, r_hz_habitable_zone |
| `calculate_quantum_gravity` | master_equation, GR_macroscopic_limit, QFT_microscopic_limit |
| `calculate_black_hole` | r_min_form_A/B/C |
| `calculate_bsd_rank_cohomology` | L_prime_E_1_uqff |
| `calculate_vds_dvp_bh26` | vds_floor, dvp_prime_anchor, bh26_92GHz_bins |
| `calculate_negative_time_dual_existence` | t_neg, CW/CCW branches |
| `calculate_si_derivations` | C_LIGHT_uqff_m_per_s, G_NEWTON_uqff_microscopic |
| `calculate_shell_orbital` | Mayer-Jensen shell + UQFF magic match |
| Primitive surfaces (5) | `calculate_resonant_adpm`, `calculate_scm`, `calculate_triadic_g`, `calculate_f_u_bi`, `calculate_f_u_bi_i` |

---

## 5. Axiom proof sets + supporting theorems = **1,867 whitepapers**

**Directory:** `whitepapers/`
**File pattern:** `PAPER_<NNNN>_<descriptive_slug>.md`
**Total files:** 1,867 (PAPER_001 through PAPER_1795, plus PAPER_S201-S205 supplements)

Each whitepaper is the **axiomatic derivation, proof, or supporting theorem** behind
one or more closures in the calculator. The closure is a one-line numerical
evaluation; the whitepaper carries:

- The axioms invoked
- The algebra / mathematical derivation
- The physics reasoning
- The references / cross-citations to other PAPER_XXXX
- The residual vs. observation
- Sometimes diagrams, tables, or worked examples

### Whitepaper organization by paper-number range

| Range | Count | Theme |
|---|---:|---|
| PAPER_001-099 | 100 | Early foundational papers (gravitational waves, GR, basics) |
| PAPER_100-499 | 401 | Core physics derivations (most of the bucket observables) |
| PAPER_500-999 | 503 | Extensions, astrophysics, AGN, BSM |
| PAPER_1000-1199 | 230 | Sector closures (Lagrangian L_YM, L_KK, L_LENR, etc.) |
| PAPER_1200-1499 | 304 | Clay Millennium proofs + paradox theorems |
| PAPER_1500-1795 | 296 | LANDMARK primitive reduction (PAPER_1521/1522), recent additions |
| PAPER_S201-S205 | 5 | Phase-Horizon supplements |

### Notable axiom-proof whitepapers

| Whitepaper | What it proves |
|---|---|
| `PAPER_001_GW170817_UQFF_Damping_Analysis.md` | GW damping derivation |
| `PAPER_646_Universal_Inertial_Operator.md` | U_i = 2.75e-7 derivation |
| `PAPER_1318_*Yang_Mills*.md` | YM mass gap = 1.736 GeV |
| `PAPER_1133_Holmlid_Rydberg_SCm_Bridge.md` | 630 eV chain derivation |
| `PAPER_1141_Rossi_ECat_Variants_Unified.md` | Star-Magic reactor closure |
| `PAPER_1156_Lambda_UQFF_Cosmology.md` | Λ at 0.003% match |
| `PAPER_1167_All_8_Lagrangian_Gaps_Closed.md` | Master synthesis |
| `PAPER_1182_*Riemann_Hypothesis*.md` | t_10000 EXACT derivation |
| `PAPER_1203_Canonical_v1.5.md` | F_U=0 master equation |
| `PAPER_1230_*Hodge*.md` | Hodge conjecture closure |
| `PAPER_1521_D_BSFG_Derivative_From_D_crit.md` | D_BSFG = 6 EXACT, primitive reduction |
| `PAPER_1522_K_MEX_Derivative_From_Phi_5_6.md` | K_MEX = 25/12 EXACT, primitive reduction |
| `PAPER_S201-S205_Phase_H*` | Final-phase gap-closure supplements |

### Browse + search whitepapers

```bash
ls whitepapers/                                    # all 1,867 files
ls whitepapers/ | wc -l                            # count
grep -l "Yang-Mills" whitepapers/*.md              # papers mentioning YM
grep -l "Holmlid" whitepapers/*.md                 # papers mentioning Holmlid
grep -l "rho_SCm" whitepapers/*.md | head          # papers using ρ_SCm

# Most-recent additions
ls -lt whitepapers/ | head -20

# By number prefix
ls whitepapers/PAPER_152*.md                       # the LANDMARK papers
```

---

## 6. C++ reference implementation

**File:** `uqff_exact_closures.cpp` (60 KB, 368 functions)

Standalone C++17 file with no dependencies beyond `<cmath>`. Each function is an
independent re-derivation of a high-value EXACT closure from UQFF.

**Purpose:** cross-language verification. If `uqff_pure_calculator.py` returns
X and `uqff_exact_closures.cpp` returns X, the closure is doubly confirmed
across two independent implementations. This is the closest thing UQFF has to
a "second opinion" within the codebase.

**Browse:**

```bash
grep "^double " uqff_exact_closures.cpp           # all 368 function signatures
```

Tier-3 K1 in PRODUCTION_ROADMAP.md is to extend this to the full 794 closures.

---

## 7. Fidelity gate — 857 runtime axiom tests

**File:** `uqff_fidelity_tests.py`
**Run:** `python uqff_fidelity_tests.py` or `uqff gate`

Breakdown:
- 354 `_exact()` regression pins — assert `|UQFF − target| / target < tol`
- 241 `check()` general assertions — structural invariants, dispatch tests
- 262 additional category-level tests across 60+ test blocks

**These are the live, executable proofs that every closure still holds its claimed value.**

If a future edit silently breaks a closure, the gate exits non-zero and the
CI matrix goes red on the next push. This is the project's anti-drift mechanism
(per CLAUDE.md Rule 8).

---

## Discovery cheat sheet

```bash
# Browse / search across ALL namespaces (after `pip install uqff`)
uqff list --all                          # all ~1,080 named closures
uqff list --all --filter <substr>        # filter
uqff search <substr>                     # search

# Look up a specific closure
uqff predict <name>                      # try all 5 namespaces

# Run the proofs
uqff gate                                # 857-test fidelity check
python uqff_fidelity_tests.py            # same

# Status summary
uqff status                              # production-readiness summary
uqff version                             # version + headline metrics

# Whitepapers (axiom proofs / supporting theorems)
ls whitepapers/                          # 1,867 PAPER_XXXX.md files
grep -l "Yang-Mills" whitepapers/*.md    # papers on a topic

# C++ reference
grep "^double " uqff_exact_closures.cpp  # 368 C++ functions

# Python source (the dispatch + closures)
grep "^def _l96_uqff_axiom_.*_closure" uqff_pure_calculator.py | wc -l
grep "^def _millennium_" uqff_pure_calculator.py
```

---

## What's NOT included in the calculator

Per CLAUDE.md Rules:
- No docstrings (Rule 3) — explanatory layer lives in whitepapers + this Atlas
- No comments inside closure functions (Rule 3)
- No SM-named constants (Rule 4 — UQFF doesn't share with SM)
- No metadata dict-keys beyond `{'value': X}` for public surfaces (Rule 5)
- No `datetime`, `__main__`, file writes, classes (Rule 6)

The calculator is intentionally a pure numerical-evaluation engine. The
**reasoning, axioms, and proofs** live in the whitepapers/ directory and the
fidelity gate.

---

## Cross-reference quick map

| You want to find... | Look here |
|---|---|
| "What does UQFF say about <observable>?" | `uqff search <observable>` |
| "Where's the derivation of <observable>?" | `grep -l "<observable>" whitepapers/*.md` |
| "Is this closure still passing?" | `uqff gate` (or check CI status on GitHub) |
| "What primitives does <closure> use?" | Read the whitepaper PAPER_XXXX referenced as primary_source |
| "How does UQFF compare to SM on <observable>?" | Bucket observable's `anchor` field is the measured value |
| "Why is <primitive> = <value>?" | `PROVENANCE_AUDIT.md` |
| "What's the parameter count?" | 9 truly-independent + 2 derivative; see `PROVENANCE_AUDIT.md` |
| "Show me all 8 Clay Millennium" | `uqff search hodge`, etc., or section 1 above |
| "Show me all LENR reactors" | `uqff predict <reactor_key>` from section 4 list above |
| "Show me a worked example" | `notebooks/` directory (00_quickstart, 01_holmlid_lenr, 02_magic_numbers, 03_cosmology) |
