# UQFF Calculator — Input Domain Documentation (Tier-1 B3)

**Last updated:** 2026-06-18
**Scope:** All 34 public `calculate_*` surfaces in `uqff_pure_calculator.py`.

This document specifies, for every public calculator surface: what `dataset` keys it consumes, valid value ranges, units, semantics, and what the return value means. All surfaces conform to Rule 5 of CLAUDE.md (return `{'value': X}` only).

---

## Universal contract

Every public surface follows this contract:

```python
result = calculate_<surface>(dataset: dict) -> {'value': X}
```

- `dataset` is always a Python `dict`. May be empty `{}` to get default-parameter behavior.
- Return is always `{'value': X}` where X may be a number, dict, or list-of-observables.
- No exceptions are raised for malformed input; instead `{'value': None}` is returned.
- All surfaces are stateless and side-effect-free (Rule 6 — no datetime, file writes, json.dump, classes).
- All numerical values use SI units unless explicitly noted.

---

## 1. Primitive computation surfaces (returns scalar)

### `calculate_resonant_adpm(dataset)`
**Returns:** ω·cos(π·t_n)·Φ_res (resonant ADPM frequency, Hz)
**Dataset keys:**
| key | type | default | range | unit |
|---|---|---|---|---|
| omega | float | 1.25e12 | (1e11, 1e14) | Hz |
| t_n | float | 0.0 | [0, 1] | dimensionless |
| phi_res | float | 0.84 | (0, 1) | dimensionless |

**Empty-dataset result:** `1.05e12` (Hz)

### `calculate_scm(dataset)`
**Returns:** SCm 26-level energy density with t_n modulation (J/m³)
**Dataset keys:**
| key | type | default | range | unit |
|---|---|---|---|---|
| rho_scm | float | 7.09e-37 | locked | J/m³ |
| s_26 | float | 1.453162 | locked | dimensionless |
| t_n | float | 0.0 | [0, 1] | dimensionless |

**Empty-dataset result:** `8.508e-37` (J/m³)

### `calculate_f_u_bi(dataset)`
**Returns:** Universal buoyancy F_UBi (dimensionless or N depending on G·M·r² inputs)
**Dataset keys:**
| key | type | default | range | unit |
|---|---|---|---|---|
| G | float | 6.674e-11 | (>0) | m³/(kg·s²) |
| M | float | 1.989e30 | (>0) | kg (solar mass default) |
| r | float | 1.496e11 | (>0) | m (1 AU default) |
| t_n | float | 0.0 | [0, 1] | dimensionless |

**Empty-dataset result:** `1.541e-15`

### `calculate_f_u_bi_i(dataset)`
**Returns:** F_U_Bi_i 4-layer master integral
Same dataset keys as `calculate_f_u_bi` plus:
| key | type | default | range | unit |
|---|---|---|---|---|
| omega_s | float | 2.5e-6 | (>0) | rad/s |

**Empty-dataset result:** `5.669e-24`

### `calculate_triadic_g(dataset)`
**Returns:** Triadic gravity: w_C·g_comp + w_R·g_res + w_B·g_buoy
**Dataset keys:** as `calculate_f_u_bi_i`. Returns dict `{triadic, g_comp, g_res, g_buoy}`.

### `calculate_vacuum_ledger(dataset)`
**Returns:** 4-term: V0 + R26 + ρ_KK + ρ_BSFG → Planck Λ (J/m³)
**Dataset keys:** typically empty (all default to canonical primitives).

---

## 2. Symbolic dispatch hub

### `calculate_analytic_closures(dataset)`
**Returns:** symbolic dispatch hub — routes to internal closures.
**Dataset keys:**
| key | type | default | range | unit |
|---|---|---|---|---|
| target | str | None | any closure name | — |

**Empty-dataset result:** `{'value': None}` (no target → no dispatch).

---

## 3. Canonical PAPER-1141 / PAPER-646 restored surfaces

### `calculate_universal_inertial_operator(dataset)`
**Returns:** dict with `U_i_dimensionless`, `U_i_product_J_per_m3`, breakdown.
**Dataset keys:**
| key | type | default | range | unit |
|---|---|---|---|---|
| omega_s | float | 2.5e-6 | (>0) | rad/s (stellar rotation) |
| t_n | float | 0.0 | [0, 1] | dimensionless |
| f_trz | float | 0.1 | locked | dimensionless |
| lambda_i | float | 1.0 | locked | dimensionless |

**Empty-dataset result:** `U_i = 2.75e-7` (Sun, t=0) — matches PAPER_646 exactly.

### `calculate_nuclear_magic(dataset)`
**Returns:** dict with `magic_numbers` (all 7), `BE_per_A_peak_Fe56`, `deuteron`, `alpha_binding`, etc.
**Dataset keys:** typically empty. Optional `Z`, `A`, `N` to focus on specific nuclide.

### `calculate_lenr(dataset)`
**Returns:** single-variant LENR closure report.
**Dataset keys:**
| key | type | default | range | unit |
|---|---|---|---|---|
| variant | str | 'holmlid' | {'holmlid','parkhomov','pons','mizuno','rossi'} | — |
| reactor_params | dict | {} | per-variant | — |

**Empty-dataset result:** holmlid KER chain (E_phonon = 8.28e-22 J ≈ 5.17 meV).

### `calculate_f_u_zero(dataset)`
**Returns:** dict with `F_UBi`, `F_UBii`, `F_UBi_plus_UBii`, `r_hz_habitable_zone` (m).
**Dataset keys:**
| key | type | default | range | unit |
|---|---|---|---|---|
| star_M_kg | float | 1.989e30 | (>0) | kg |
| Z | float | 1.0 | (>0) | atomic number |
| t_n | float | 0.0 | [0, 1] | dimensionless |
| r_test | float | 1.496e11 | (>0) | m (initial guess) |

---

## 4. Phase-2 canonical surfaces

### `calculate_ua_layers(dataset)`
**Returns:** dict `{UA_prime, UA_double, UA_triple, UA_quad, F_DPM}` (all J/m³).
**Dataset keys:** typically empty.

### `calculate_dpm_grinding(dataset)`
**Returns:** 5-step DPM grinding sequence (SCm → UA → UA''''').
**Dataset keys:** typically empty.

### `calculate_caduceus(dataset)`
**Returns:** 26 Caduceus pinch points + π decimal encoding.
**Dataset keys:** typically empty.

### `calculate_shell_orbital(dataset)`
**Returns:** Mayer-Jensen shell occupancy at given Z + UQFF magic-number match.
**Dataset keys:**
| key | type | default | range | unit |
|---|---|---|---|---|
| Z | int | 2 | [1, 126] | atomic number |

### `calculate_lenr_full(dataset)`
**Returns:** Full LENR report — 8 reactors + Widom-Larsen + Kozima + meson cascade.
**Dataset keys:** typically empty.

---

## 5. BUCKET 0 loop-closure surfaces

### `calculate_negative_time_dual_existence(dataset)`
**Returns:** t_neg < 0 + CW/CCW dual branches.
**Dataset keys:** typically empty.

### `calculate_si_derivations(dataset)`
**Returns:** dict with `C_LIGHT_uqff_m_per_s` and `G_NEWTON_uqff_microscopic`, parameter-free.
**Dataset keys:** none.

### `calculate_quantum_gravity(dataset)`
**Returns:** 26D quantum-gravity unification (master equation, GR+QFT limits, mass-gap, no-singularity).
**Dataset keys:** typically empty.

### `calculate_black_hole(dataset)`
**Returns:** dict with `r_min_form_A_eigenvalue_m` and `r_min_form_B_triad_m` per mass.
**Dataset keys:**
| key | type | default | range | unit |
|---|---|---|---|---|
| M_bh_kg | float | 1.989e30 | (>0) | kg (solar default) |

### `calculate_vds_dvp_bh26(dataset)`
**Returns:** VDS (P/3 floor) + DVP (prime=113) + BH26 (92 GHz × 26 bins).
**Dataset keys:** typically empty.

### `calculate_bsd_rank_cohomology(dataset)`
**Returns:** dict with `L_prime_E_1_uqff` (0.30598, 0.005%) vs anchor.
**Dataset keys:** typically empty.

---

## 6. BUCKET A Millennium dispatch

### `calculate_paradox(dataset)`
**Returns:** routes paradox closure via PARADOX_TO_CLOSURE dispatch.
**Dataset keys:**
| key | type | default | range | unit |
|---|---|---|---|---|
| paradox | str | None | any of 802 paradox names (lowercase, see PARADOX_TO_CLOSURE) | — |

**Critical:** the `paradox` value is normalized via `.lower().strip().replace("-","_").replace(" ","_")`. Use **lowercase** dispatch keys. See CLAUDE.md case-sensitivity note.

**Empty-dataset result:** `{'value': {'inventory': {'total_paradoxes': 802, ...}}}` — describes the dispatcher state, not a closure.

---

## 7. BUCKET C-K observable suites

These all return `{'value': {'observables': [...]}}` where each observable is `{'observable': label, 'uqff_derived': val, 'anchor': anchor, 'residual_pct': r}`.

### `calculate_cosmology(dataset)` — 56 cosmological observables
**Dataset keys:** typically empty. Optional `redshift_z`, `H0_anchor_km_s_mpc` for hyperscale variants.

### `calculate_particle_physics(dataset)` — 48 particle-physics observables
**Dataset keys:** typically empty.

### `calculate_gw_events(dataset)` — 22 gravitational-wave events
**Dataset keys:** typically empty. Optional `event_name` to focus on one event.

### `calculate_agn_jet(dataset)` — 22 AGN/jet observables
**Dataset keys:** typically empty.

### `calculate_astrophysics(dataset)` — 36 astrophysics observables
**Dataset keys:** typically empty.

### `calculate_high_energy_astro(dataset)` — 10 TeV/PeV observables
**Dataset keys:** typically empty.

### `calculate_qgp(dataset)` — 9 QGP observables
**Dataset keys:** typically empty.

### `calculate_higgs_precision(dataset)` — 13 Higgs observables
**Dataset keys:** typically empty.

### `calculate_bsm_constraints(dataset)` — 17 BSM observables
**Dataset keys:** typically empty.

---

## 8. Tier-1 utility surfaces

### `calculate_status_report(dataset)`
**Returns:** dict with `summary` (counts, classifications, milestones).
**Dataset keys:**
| key | type | default | range | unit |
|---|---|---|---|---|
| include_bucket_breakdown | bool | False | True/False | — |

**Empty-dataset result:** comprehensive summary including:
- total_closures: 794
- uncertainty_classes_A2_TIER1_production_readiness (5-band classification)
- with_full_schema: 263
- legacy_freeform: 530
- unique_paper_sources: 244
- cosmic_milestones (5 booleans)
- truly_independent_primitives: 9
- derivative_primitives: 2

### `calculate_whitepaper(dataset)`
**Returns:** whitepaper catalog summary.
**Dataset keys:** typically empty.

---

## Common pitfalls

1. **Mixed-case dispatch keys silently fail** for `calculate_paradox`. The dispatcher lowercases input but does NOT raise on miss. Always use lowercase keys like `"hubble_tension"`, `"strong_cp"`, `"vacuum_catastrophe"` — never `"Hubble_Tension"`.

2. **Empty dataset returns the canonical-primitive default**, not None. To get a "no data" signal, set explicit `target=None` or `paradox=None`.

3. **Time `t_n` is dimensionless in [0,1]**, not in seconds. Multiply by canonical period if you need a physical-time interpretation.

4. **All masses in kg, distances in m, frequencies in Hz**. There is no internal unit-conversion layer.

5. **Returned `'anchor'` value is the literal measured observable**, not a SM-theory prediction. Per Rule 4, the calculator never compares to SM-derived anchors; only to direct experimental measurements.

6. **Return-value dicts may be deeply nested**. Use `.get('value', {})` chains, not direct attribute access.

7. **The 530 legacy_freeform closures only have `{'value': X}` — no anchor/residual.** Schema-tagged closures (263) carry the full `(target, residual_pct, status_tier, description, primary_source)` tuple.

---

## Reference example

```python
import uqff_pure_calculator as u

# Get the Holmlid KER prediction
result = u.calculate_lenr({'variant': 'holmlid'})
ker_eV = result['value']['ker_chain']['E_phonon_eV']
# (Should be ~5.17 meV per phonon, scaled by S_26 and Φ_res to 630 eV total)

# Compute habitable zone for a 2-solar-mass star
result = u.calculate_f_u_zero({'star_M_kg': 2 * 1.989e30, 'Z': 1.0, 't_n': 0.0})
r_hz_m = result['value']['r_hz_habitable_zone']
# Returns r_hz in meters

# Look up cosmological observable suite
result = u.calculate_cosmology({})
for obs in result['value']['observables']:
    print(f"{obs['observable']}: UQFF={obs['uqff_derived']:.4e}, anchor={obs['anchor']:.4e}, residual={obs['residual_pct']:.3f}%")

# Look up a paradox closure (use LOWERCASE dispatch key!)
result = u.calculate_paradox({'paradox': 'hubble_tension'})
value = result['value']

# Get the comprehensive status report
result = u.calculate_status_report({})
print(result['value']['summary'])
```

---

**Bottom line for production reporting:**

> All 34 public `calculate_*` surfaces follow a stateless `dict -> {'value': X}` contract. Most accept an empty `dataset` and return canonical-primitive defaults. The 9 bucket-observable surfaces (cosmology, particle physics, GW, AGN, astrophysics, high-energy astro, QGP, Higgs, BSM) return structured observable lists with anchor + residual. The single `calculate_paradox` dispatch surface routes 802 paradox closures via a lowercase-normalized key. No surface raises exceptions — malformed input returns `{'value': None}`. Tier-2 packaging will expose this API surface via `uqff.predict(...)`.
