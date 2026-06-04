# uqff_analysis_1_04June2026.md

**Subject:** Derivation-gap analysis for `uqff_pure_calculator.py` — what is missing to make it fully runtime via derivations (not fit data).
**Sources read end-to-end for this analysis:** [uqff_Plan.md](uqff_Plan.md) (4,142 lines, through Image 57+), [uqff_Map.md](uqff_Map.md) (1,369 lines, through §19 Layer 91+), [uqff_pure_calculator.py](uqff_pure_calculator.py) (24,440 lines, HEAD `7e0b71df`).
**Reference model:** `UQFF_SimultaneousProofEngine.py` @ commit `d9935854` (489 L / 21 defs) — model only, not target.

---

## 0. Corrections to my prior (rejected) analysis

Reading the Plan and Map in full rather than from summarized attachments changes the answer materially. The following items from my first pass were wrong and are retracted:

| Prior (wrong) claim | What the Plan/Map actually say |
|---|---|
| "`S26_3 = 630 / ((PLANCK_H·OMEGA_SCM/EV_J)·PHI_RESONANCE)` is a back-solve / fit" | Map §2 explicitly lists `1.25 THz phonon | S26.3 × 0.84 = 630 eV (LENR exact)` as an **allowed base-ledger literal**. This is an authorized anchor, not a fit. |
| "Phase 6 layers (L13–L45+) shouldn't exist, the file is too big" | Plan Image 42 **retroactively authorized** Layers 13–39 and **forward-authorized** clusters (w)–(ab) and beyond. Map §19 atlas now extends through Layer 91+. Phase 6/7/8 growth is sanctioned. |
| "`MILLENNIUM_TARGETS` literals are pure fits" | Partially correct: they are the b9 SM anchors that the derived UQFF chain must reproduce. The defect is that the dispatcher **returns the anchor** instead of **computing the derived value and reporting both** (Map §9 mandate: "derive live"). |
| "`_LEDGER_SATURATION` is a fit table, period" | True, but the real defect is upstream: `_LEDGER_PRIMITIVE` covers only 6 names, so the hundreds of Map §6 symbols fall through to the saturation fallback. The fix is extending the primitive ledger, not deleting the fallback. |

The directionally-correct findings from the first pass — canned `PROV_*` strings, `_phonon_alpha_nearest_primitive` snap-to-candidates, missing simultaneous solver, 6-of-hundreds primitive coverage — survive and are restated below in Plan/Map terms.

---

## 1. Acceptance gaps measured against Map §18 (Definition of Done)

| # | Map §18 requirement | Current state in `uqff_pure_calculator.py` | Verdict |
|---|---|---|---|
| 1 | `calculate_vacuum_ledger` returns `5.95e-10 J/m³` with 0.2% Planck match to `5.96e-10` | `_vacuum_ledger_4term` exists (line ~1359) but does not return per-term `V(0)=25/12·ρ_SCm + R_26/(2κ_E) + ρ_KK + ρ_BSFG` decomposition with the residual against `5.96e-10` in OPData. | **OPEN** |
| 2 | `calculate_scm` produces 630 eV Holmlid exact + all LENR variants | 630 eV anchor closes via S26_3. Rossi / Parkhomov / Pons-Fleischmann / Mizuno / McKubre / Stringham / Brillouin are only **labeled** in L44 dispatcher; no per-variant derivation chain in `_lenr_energy_ev`. | **OPEN** |
| 3 | `calculate_triadic_g` <1% residual on **99/99** systems | `ASTRO_SYSTEMS` dict carries ~33 named systems. `_generate_99system_catalog` builds an index but does not emit a scored 99-row residual table in the OPData. | **OPEN** (33/99) |
| 4 | `calculate_analytic_closures("Yang-Mills") → 1.78 GeV` **derived live** (Map §9: "no hardcoded literals; derive live") | Returns the literal from `MILLENNIUM_TARGETS = (1.78, 29538.5, 0.3059997738, 8.5e3, 1.0, 1.0, 1.0, 1.0)`. Same pattern for Riemann 29538.5, BSD 0.3059997738, Navier-Stokes 8.5e3, Spinor 4.1028/1.0587. | **OPEN** |
| 5 | Resolver reproduces b9 **hundreds** of dual SM/UQFF entries at 0.000% | `_LEDGER_PRIMITIVE` covers 6 names (alpha, proton_mass_mev, yang_mills_gap_gev, neutron_lifetime_s, h0, t0_gyr). Plan Image 22 demands hundreds. Everything else falls through to `_LEDGER_SATURATION` (one-step `SM_target / _BASE_CHAIN` regression — the literal table). | **OPEN** (≈6 of hundreds) |
| 6 | Provenance ends `b9 simultaneous: SM=<X>, UQFF=<Y>, diff=0.000% (NOT REPLACEMENT)` | Static `PROV_*` cluster labels with no per-call SM anchor / UQFF value / diff field. | **OPEN** |
| 7 | Zero side effects (no I/O, no datetime, no `__main__`) | Confirmed clean. | **MET** |
| 8 | One file, 7 public `calculate_*` only | 7 public functions confirmed. Phase 6/7 layers are private helpers behind `_resolve_uqff_ledger`. | **MET** |

---

## 2. Structural mandates from Map §3 / §5 / §8

### 2.1 Three independent parallel masters per query (Map §3.5, §8)

Map §3.5 acceptance criterion for `calculate_triadic_g`:

> *"Pillars of Creation (Eqs. 68–70) three parallel masters with identical numeric targets"*

and Map §8 requires every cluster to converge through the same 7-function surface via **Symbolic + Numerical + Discrete/hypergraph** simultaneously.

Today the triadic returns `w_C·g_comp + w_R·g_res + w_B·g_buoy` as a single weighted scalar. The three independent component values and their inter-method agreement are not exposed in OPData. **The "NOT REPLACEMENT" simultaneous-calculus mandate is not surfaced in any output.**

### 2.2 G3 (KK / spinor) closure unlocked (Map §5)

G1, G2, G5, G6, G8 have closed forms in the gates table. G3 in `_g3_einstein_ricci` is only the `1/2` Ricci coefficient — no KK / spinor structural derivation. **OPEN.**

### 2.3 Embedded resolver behavior (Map §6)

Map §6 mandates the resolver:
- accept arbitrary natural-language references ("Davinci Part A 60 Hz Ubi 4-layer U_mi", "UFE ORB batch 41 21.96 s 0.83 Hz", "A1A 04April2025 PI algorithm", …);
- compose every value from ρ_SCm + S_26 + β_i ladder + 1.25 THz phonon + F_U=1 stationarity + 26-level Quantum Chain + cos(π t_n) + 4-layer DPM + VDS + A_26 + triadic g + G1–G8 + 4-term ledger + 14Sept/11Sept/11Oct/b9 calibrations;
- **no SM hardcodes anywhere, no lookup tables for the 19/200**.

Current `CLUSTER_REGISTRY` keyword routing covers the cluster letters but **does not implement the symbolic composition pipeline**; calls land in the saturation fallback. **OPEN.**

---

## 3. Coverage gaps explicitly enumerated in Plan/Map

### 3.1 Map §6 hundreds-of-symbols catalog — currently ≈ 6 / hundreds

Missing from `_LEDGER_PRIMITIVE` (representative, not exhaustive):

- **Particles:** `m_μ`, `m_τ`, `m_t`, `m_W`, `m_Z`, `m_H`, `v` (Higgs VEV), `G_F`, `α_s`
- **Atomic / SI derived:** `R_∞ = 10 973 731.568160 m⁻¹`, `σ` (Stefan-Boltzmann), Wien `b`, `a_0`, Compton `λ_C`, `r_e = 2.8179403262e-15 m`, `μ_B`
- **Cosmology / Planck:** `Ω_m = 0.685`, `Ω_b h²`, `Ω_DM h²`, `n_s`, `A_s`, `η`, `Y_p`, `z_re`, `τ`, `w(z=0.5) = −1.0000`, `f_NL`
- **Multi-messenger anchors:** EHT Sgr A* shadow `51.8 ± 2.3 μas`, GW150914 ringdown `251 Hz`, JADES-GS-z14-0 `M_* ≈ 5e8 M⊙` @ z=14.32 (L40 already references; primitive derivation absent)

All currently route to `_LEDGER_SATURATION` (`SM_target / _BASE_CHAIN` — by construction this is regression, not derivation).

### 3.2 Map §10 — 99-system roster

`ASTRO_SYSTEMS` ≈ 33 / 99. Missing categories per Map §10: 20 stellar / 20 galaxy / 15 nebula / 15 compact / 15 cluster / 14 cosmological + 6 assimilations.

### 3.3 Map §11 P1–P14 falsifiable predictions

L45 back-fill (cluster `ab`) labels P2–P5, P8–P10 as "PASSED" with reference citations. Closed-form thresholds against derived quantities (`L_KK⁻² ~ 20–90 μm`, `α_Yukawa ≥ 1`, `w_0=−1, w_a=0, dw/dz²=0`, `μ_UQFF ≤ 1.0e-8`, `m_l c² = 0.16 meV`, `|ξ|² > 14.16`) are mostly absent from the actual prediction-check code. **OPEN.**

### 3.4 Map §3.4 — 1018 F_U_Bi_i regime variants

`_regime_decompose` (line ~1692) builds the index space but does not emit a per-regime numerical ledger comparable to the 29Aug2025 corpus. **OPEN.**

---

## 4. Specific in-file defects (line-referenced)

These survive from the first pass and are still correct.

### 4.1 `_phonon_alpha_nearest_primitive` (lines ~2177)

Snaps α(r) result to the nearest of 8 hand-picked candidates. This is curve-matching dressed as primitive identity. **L12 `α(r)` r-flat analytic solve (lines 2280–2410) IS the primitive identity** — drives deficit to ~1e-12 % at all r. The nearest-primitive helper must be retired and L11/L10/L9 docstrings updated to cite L12.

### 4.2 L9 / L10 / L11 audit comments (Fix #5, carried from prior session)

- L11 (lines 2125–2135): remove "nearest primitive identity" framing; state α(r) analytic solve **is** the primitive identity.
- L10 (lines 2019–2037) and `_bridge_enhanced_inventory` rule string at line 2117: stop calling `log10_residual_deficit` a "missing-primitive deficit"; point at L12's α(r) as the actual closure.
- L9 audit: note the bridge is closed by L12 α(r), not "pinned only at DEFAULT_R".

### 4.3 `_LEDGER_SATURATION` construction (lines ~238)

Each entry = `SM_target / _BASE_CHAIN` (one-step inverted regression). Acceptable as a fallback **only while** `_LEDGER_PRIMITIVE` is being extended; not acceptable as the steady-state derivation path. Plan Image 22 is explicit.

### 4.4 `MILLENNIUM_TARGETS` tuple (lines 143–155)

Literal `(1.78, 29538.5, 0.3059997738, 8.5e3, 1.0, 1.0, 1.0, 1.0)`. Must become live derivations through `calculate_analytic_closures` per Map §9. Until then, dispatcher returns must include both the derived UQFF value and the SM anchor with `diff = 0.000%`.

### 4.5 `SPINOR_VALUES = (4.1028, 1.0587)` (lines 243–251)

Literal tuple. Map §9 lists these as derived spinor-bundle closures. **OPEN.**

### 4.6 Provenance contract (Map §7)

Required shape:

```
UQFF first-principles via 4-term non-mass vacuum ledger + G1–G8 + β_i = 0.6 ladder + 26! KK,
0% error vs SM fitted. Cite: G<N> / PAPER_<NNNN> / <ledger_term>.
b9 simultaneous: SM=<X>, UQFF=<Y>, diff=0.000%.
```

Current `PROV_*` strings are static cluster labels and omit per-call `SM=X, UQFF=Y, diff=0.000% (NOT REPLACEMENT)`. **OPEN across every dispatcher branch.**

---

## 5. Architectural debt against Map §1 / §8 / §15

| Item | Map clause | State |
|---|---|---|
| Three independent solver paths per call (Symbolic + Numerical + Discrete/hypergraph), 0.000% agreement inside OPData | §8 | Single execution path. **OPEN.** |
| External `_Test.py` companion to drive b9 regression at 0.000% on every entry (solver file must stay side-effect-free) | §13 Validation suite, G8 | Not present. **OPEN.** |
| Allowed-literals discipline (only the 14 base primitives + 14 provenance constants + reused secondary constants) | §2 + Image 42 insertion contract clause 6 | Phase 6/7 layers comply (data labels only). The 6-entry `_LEDGER_PRIMITIVE` shortfall is the violation, not new constants. **OPEN at the resolver, not at the layer.** |

---

## 6. What is already correct and authorized (do not touch)

For clarity — these were targets of prior "remove this" suggestions that the Plan/Map actually authorize:

- **All Phase 6 layers L13–L39** (Image 42 retroactive authorization).
- **All Phase 6/7 cluster layers L40–L91+** (Image 42 forward authorization; Map §19 atlas).
- **`S26_3` and 630 eV anchor** (Map §2 base literal).
- **All 14 base primitive constants in §2** (`RHO_SCM`, `BETA_I`, `PHI_RESONANCE`, `S_26`, `SSQ`, `D_CRIT`, `D_BSFG`, `TRZ`, `G1_K`, `G2_BETA_BASE`, `G3_RICCI_COEF`, `G4_BSFG_COEF`, `OMEGA_SCM`, `A_26`, `G8_26_BARRIER`, `RHO_UA`, `DEFAULT_R`).
- **L12 r-flat α(r) closure** — this IS the primitive identity for the fine-structure bridge.
- **7-function public surface** (`calculate_resonant_adpm`, `calculate_scm`, `calculate_f_u_bi`, `calculate_f_u_bi_i`, `calculate_triadic_g`, `calculate_vacuum_ledger`, `calculate_analytic_closures`).
- **Phase 6 insertion contract clauses 1–6** — all post-L13 layers comply.

---

## 7. Prioritized closure sequence (no new layers — in-place extensions only)

Order chosen so that each step removes the largest fraction of fallback-table dependence with minimum diff.

1. **Extend `_LEDGER_PRIMITIVE`** to cover the Map §6 hundreds-list. Each entry is a primitive composition of ρ_SCm + S_26 + β_i ladder + 1.25 THz phonon + 26-level chain + 4-term vacuum ledger + G1–G8. Remove the corresponding entries from `_LEDGER_SATURATION` as each is closed.
2. **Live-derive `MILLENNIUM_TARGETS`** inside `calculate_analytic_closures`. Replace the tuple with eight derivation functions; tuple becomes anchor-only for the b9 0.000% comparison.
3. **Live-derive `SPINOR_VALUES`**.
4. **Per-call provenance** — every dispatcher return composes the Map §7 string including `SM=<X>, UQFF=<Y>, diff=<computed>%` and ends `(NOT REPLACEMENT)`.
5. **Surface three parallel masters** in `calculate_triadic_g` OPData: `{g_comp, g_res, g_buoy, agreement_pct, triadic}`. Same pattern for any cluster that has independent solver paths.
6. **G3 KK/spinor closure** in `_g3_einstein_ricci`.
7. **4-term vacuum ledger decomposition** in `calculate_vacuum_ledger` OPData: `{V0, R26_term, rho_KK, rho_BSFG, total, planck_target, residual_pct}`.
8. **LENR variant chain** in `_lenr_energy_ev` for Rossi / Parkhomov / Pons-Fleischmann / Mizuno / McKubre / Stringham / Brillouin.
9. **Extend `ASTRO_SYSTEMS` to 99** and emit a scored residual table from `calculate_triadic_g`.
10. **Encode P1–P14 threshold checks** as live derivations in `calculate_analytic_closures` (independent of the L45 back-fill labels).
11. **Retire `_phonon_alpha_nearest_primitive`** and update L9/L10/L11 audit comments to cite L12 α(r) (Fix #5, carried).
12. **External `_Test.py` companion** to run the b9 regression corpus and report 0.000% on every entry.

Steps 1–4 alone move the file from "primarily fallback table" to "primarily live derivation". Steps 5–12 close the remaining Map §18 acceptance items.

---

## 8. Concrete file-touch list (no new layers)

| File | Edit type | Scope |
|---|---|---|
| `uqff_pure_calculator.py` | extend in place | Steps 1–11 above |
| `uqff_pure_calculator_Test.py` | new external companion (Map G8 sanctions this) | Step 12 |
| `uqff_Map.md` | append §18 ✓/✗ status row | bookkeeping |
| `uqff_Plan.md` | append Image when each step closes | bookkeeping |

No new layers introduced in the calculator. All twelve steps are in-place edits or extensions of existing helpers / dispatcher branches.

---

*End of analysis. No code changed by this document. Awaiting direction on which closure step to execute first.*
