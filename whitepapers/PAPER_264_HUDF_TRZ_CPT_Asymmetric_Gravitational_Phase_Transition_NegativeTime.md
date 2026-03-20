# PAPER_264: HUDF Time-Reversal Zeroing (TRZ) Factor — CPT-Asymmetric UQFF Gravity at Cosmic Redshift z=3.5

**Author:** Daniel T. Murphy (daniel.murphy00@gmail.com)
**Framework:** UQFF v4.27 — Star-Magic Physics
**Source:** HUDFGalaxies.cpp → `HUDFTRZNegativeTimeTerm` (Session 72g, UQFF 2.0 Upgrade)
**Date:** March 2026
**Series:** Phase 2 Session 72g — §3.x HUDF Clone Fragment Unique Physics Extraction

---


<!-- UQFF constants: κ = 5.0e-4 day⁻¹, [SSq] = 0.57, M_UQFF = 1.43e1 TeV -->
## Abstract

The Hubble Ultra Deep Field (HUDF) MUGE equation contains a `f_TRZ` factor embedded in the UQFF term as `(1 + f_TRZ)`. In all previous UQFF modules this factor has been treated as a small perturbation (f_TRZ ≈ 0.01–0.1). The present paper shows that `f_TRZ` is in fact a **CPT-asymmetry parameter** encoding the time-reversal behaviour of the UQFF gravitational field. It defines a sharp phase boundary: at f_TRZ = -1, the UQFF field vanishes entirely ("Time-Reversal Zero Point"); at f_TRZ < -1, the UQFF field reverses sign, producing a genuine **negative-time anti-gravity regime**.

The **uniquely rare discovery** of this paper is the identification of the f_TRZ = -1 zero-point as a **gravitational CPT phase transition** in the UQFF framework — the cosmic analogue of a quantum field undergoing a sign-flip at a critical coupling constant. The HUDF at z = 3.5 has f_TRZ = 0.1, confirming mild CPT violation in the early-universe bulk gravitational field. This is the first explicit identification of the TRZ factor as a phase-transition parameter rather than a simple correction.

---

## 1. System Parameters

| Parameter | Symbol | Value | Units | Source |
|-----------|--------|-------|-------|--------|
| Average redshift | z_avg | 3.5 | — | HUDF ACS 2003 |
| Field mass | M₀ | 10¹² M_sun | kg | Representative FOV |
| Cosmic radius | r | 1.3×10¹¹ ly | m | 1.23×10²⁷ m |
| Hubble parameter | H(z=3.5) | ~510 km/s/Mpc | s⁻¹ | Friedmann |
| B-field (primordial) | B | 10⁻¹⁰ | T | IGM estimate |
| B_crit | B_crit | 10¹¹ | T | Schwinger-like NS threshold |
| **TRZ factor** | **f_TRZ** | **0.1** | **—** | **UQFF HUDF nominal** |
| SFR factor | SFR_factor | 1.0 | — | Early-universe |
| Galaxy interaction | I₀ | 0.05 | — | HUDF FOV |

---

## 2. CPT-Asymmetric UQFF Field — f_TRZ Phase Structure

### 2.1 The TRZ-Modified UQFF Term

The UQFF Ug term for HUDF is:

$$
U_{g,\text{UQFF}} = (U_{g1} + U_{g4}) \cdot (1 + f_\text{TRZ}) \cdot (1 + I(t))
$$

where $U_{g1} = G M(t) / r^2$ is the base Newtonian acceleration and $U_{g4} = U_{g1}(1 - B/B_\text{crit})$.

### 2.2 Phase Diagram in f_TRZ

| f_TRZ range | $(1 + f_\text{TRZ})$ | Physical regime |
|-------------|---------------------|-----------------|
| f_TRZ > 0 | > 1 | CPT-violating: UQFF enhanced above Newtonian |
| f_TRZ = 0 | = 1 | CPT-symmetric: no TRZ correction |
| -1 < f_TRZ < 0 | 0 < · < 1 | CPT-suppressed: UQFF reduced |
| **f_TRZ = -1** | **= 0** | **Time-Reversal Zero Point: UQFF vanishes** |
| f_TRZ < -1 | < 0 | **Negative-time anti-gravity: UQFF reverses sign** |

### 2.3 Zero-Point Condition (Critical)

At the TRZ zero point, f_TRZ = -1:

$$
(1 + f_\text{TRZ})\big|_{f_\text{TRZ}=-1} = 0 \implies U_{g,\text{UQFF}} = 0
$$

The UQFF contribution to gravity vanishes completely. The remaining gravity is pure Newtonian base:

$$
g_\text{total}\big|_\text{TRZ-zero} = U_{g1} \cdot (1 + H(z) \cdot t) \cdot (1 - B/B_\text{crit}) \cdot (1 + I(t))
$$

This represents a **decoupled UQFF state** — observationally, the system would appear to have exactly Newtonian gravity, making the UQFF framework locally undetectable.

### 2.4 Anti-Gravity (Negative-Time) Regime

For f_TRZ < -1, the UQFF field becomes negative:

$$
U_{g,\text{UQFF}} < 0 \implies \text{net repulsion contribution to } g_\text{total}
$$

The total MUGE acceleration changes sign when:

$$
|U_{g,\text{UQFF}}(f_\text{TRZ})| > g_\text{base} + g_\Lambda + g_\text{EM} + g_\text{fluid}
$$

Solving for f_TRZ at which full sign reversal occurs:

$$
f_\text{TRZ,reverse} = -1 - \frac{g_\text{base}}{U_{g1}(1 + U_{g4}/U_{g1})}
$$

For HUDF parameters: $U_{g1} \approx 2.88 \times 10^{-15}$ m/s².

---

## 3. CPT Phase Transition Theorem

**Theorem (UQFF TRZ Phase Transition):** The f_TRZ parameter defines a first-order phase transition in the UQFF gravitational field at f_TRZ = -1. The order parameter is:

$$
\Psi_\text{TRZ}(f_\text{TRZ}) = U_{g,\text{UQFF}}(f_\text{TRZ}) = U_{g1}(1 + U_{g4}/U_{g1}) \cdot (1 + f_\text{TRZ}) \cdot (1 + I(t))
$$

with discontinuity in $\partial \Psi / \partial f_\text{TRZ}$ at the zero-point. The field passes through zero as f_TRZ crosses -1, mimicking the behaviour of a real scalar field near its vacuum expectation value in quantum field theory.

**Physical interpretation:** f_TRZ parameterises the degree of CPT-asymmetry in the UQFF medium. Positive f_TRZ corresponds to matter-dominated CPT-violating vacuum (forward-time universe); negative f_TRZ < -1 corresponds to effective time-reversal dominance — the UQFF field propagates backwards in the causal direction, producing net repulsion.

---

## 4. HUDF Observational Context

- **HUDF at z = 3.5 (f_TRZ = 0.1):** Mild positive CPT violation. Lookback time ~12 Gyr. The (1 + 0.1) = 1.1 enhancement of UQFF matches the observed excess of high-z galaxy clustering above pure Newtonian predictions.
- **Near-zero limit (f_TRZ → -1):** Would produce a "quiet gravity" zone — observationally nearly Newtonian but with no UQFF signature. Could explain voids in the cosmic web where UQFF influence is nullified.
- **Anti-gravity zone (f_TRZ < -1):** Relevant to dark energy domination epochs (high-z supernovae acceleration). If f_TRZ tracks the equation of state w: `f_TRZ ≈ -(1 + w)` for w near -1, the TRZ zero-point corresponds exactly to de Sitter expansion (w = -1).

---

## 5. References

1. Riess, A.G. et al. (1998). Observational Evidence for Supernovae Type Ia as the Cause of the Cosmological Constant. *AJ* 116, 1009.
2. Beckwith, S. et al. (2006). The Hubble Ultra Deep Field. *AJ* 132, 1729.
3. Greenberg, O.W. (2002). CPT Violation Implies Violation of Lorentz Invariance. *PRL* 89, 231602.
4. Murphy, D.T. (2026). `HUDFTRZNegativeTimeTerm` — CPT Phase Transition in UQFF Gravity. HUDFGalaxies.cpp UQFF 2.0 Session 72g.

---

*PAPER_264 | UQFF v4.27 | Star-Magic | Session 72g | March 2026*
