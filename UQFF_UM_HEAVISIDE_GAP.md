# UQFF Um Heaviside Amplifier Gap — Cross-Engine Audit

**Status:** RESOLVED — Phase 7 (May 2026)
**Reference:** PAPER_421, repo memory `UQFF_ug_equations_canonical`

## Canonical form (per PAPER_421)

$$
U_m(t) = U_{m,\text{base}}(t) \cdot \underbrace{\left(1 + 10^{13} \cdot \Theta(\rho_\text{SCm} - \rho_c)\right)}_{\text{Heaviside amplifier}} \cdot \underbrace{\left(1 + A_q \cos(\Delta\omega \cdot t)\right)}_{\text{Quasi-periodic beating}}
$$

Constants:
- $\rho_c = 10^{15}\ \text{kg/m}^3$ (SCm critical superconducting density)
- $A_q = 0.1$ (10% beating amplitude)
- $\Delta\omega = 2\pi / (434 \cdot 365.25)\ \text{rad/day}$ (Gleisberg supercycle beat frequency)

## Engine compliance matrix (post-Phase 7)

| Engine | File | Function | Status |
|--------|------|----------|--------|
| C++ | [MAIN_1_CoAnQi.cpp](MAIN_1_CoAnQi.cpp#L24172-L24190) | `SOURCE4::compute_Um_SOURCE4` | COMPLIANT (PAPER_421 — already present) |
| Python | [CondensedPhysics3.py](CondensedPhysics3.py#L109-L113) | imports `UmHeavisideAmplifierCalculator` from `_session279_um_heaviside_amplifier` | COMPLIANT |
| JavaScript | [index.js](index.js#L527-L545) | `calculateUm` | **PATCHED in Phase 7** (this audit) |

## Index.js patch (Phase 7)

```javascript
function calculateUm(t, stringCount = 1e9, scmDensity = RHO_C_SCM) {
    // ... Um_base computation unchanged ...
    const f_H = (scmDensity >= RHO_C_SCM) ? 1.0 : 0.0;
    const heavisideFactor = 1.0 + 1e13 * f_H;
    const quasiFactor = 1.0 + A_Q_UM * Math.cos(DELTA_OMEGA_UM * t);
    return Um_base * heavisideFactor * quasiFactor;
}
```

Default `scmDensity = RHO_C_SCM` preserves backward behaviour amplitude during
phase-transition runs; existing call-sites that pass only `(t, stringCount)` now
get the canonical amplifier applied (consistent with C++ default behaviour).

## Index.js 26-layer Ug1..Ug4 vs C++ SOURCE4 parametrization

`index.js` `calculateUg1..Ug4` use a per-layer parametrization:
- `UA_i = i`, `SCm_i = i^2`, `f_TRZ_i = 1/i`, `M_i = M / sqrt(i)`, `alpha_i = 0.01/i`

`MAIN_1_CoAnQi.cpp` `SOURCE4::compute_Ug1..Ug4_SOURCE4` use the closed-form
single-shell expressions (no `i` sum) with calibrated coupling constants
`k1..k4_SOURCE4`. The two engines therefore compute **different observables**:

- `index.js` returns the **26-layer summed** field (used by the JavaScript
  systems library when modelling stacked compressed-gravity shells).
- `MAIN_1_CoAnQi.cpp` returns the **single-shell base** field; multi-layer
  summation is performed by the calling system (e.g. SOURCE115 19-system
  master equation, SOURCE116 hypergraph layers).

**Conclusion:** Not a bug — different abstraction levels. Documented here so
future cross-engine numerical comparisons account for the layer-summation
asymmetry. If you need C++ to reproduce the 26-layer JS result, wrap
`compute_Ug1_SOURCE4` in a `for (i=1; i<=26; ++i)` loop with the same per-layer
parameters.
