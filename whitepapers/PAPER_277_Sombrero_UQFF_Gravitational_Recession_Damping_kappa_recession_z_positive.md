# PAPER_277 — UQFF Gravitational Recession Damping Factor κ_recession for Positive Redshift

**Author:** Daniel T. Murphy
**Module:** SOMBRERO_UQFF_MODULE.cpp (UQFF 2.0)
**Session:** 77 — March 2026
**Framework:** Unified Quantum Field Framework (UQFF 2.0)
**Status:** Complete — embedded in SOMBRERO_UQFF_MODULE.cpp
**Whitepaper Series:** 277/1000
**DOI (Provisional):** UQFF-2026-277-RECESSION

---

## Abstract

We introduce the UQFF Gravitational Recession Damping Factor κ_recession = 1/(1+z) for galaxies receding from the observer (positive cosmological redshift z > 0). Applied to the Sombrero Galaxy (M104, z = +0.0063), this factor yields κ_recession = 0.99374, attenuating the total UQFF gravitational output by 0.63% relative to the rest-frame value. This result is the precise complement of PAPER_273's blueshift amplifier (Andromeda, z = −0.001, κ > 1), and together the two whitepapers establish the **Universal UQFF Bidirectional Redshift Law**: the single analytic function κ(z) = 1/(1+z) applies for all z, unifying recession and approach within one UQFF correction framework.

---

## 1. Motivation and Context

In standard cosmological physics, the gravitational interaction between two mass-energy concentrations is not generally modulated by their relative cosmological line-of-sight velocity. However, within the UQFF (Unified Quantum Field Framework), which models gravity as a buoyancy-mediated emergent phenomenon in the quantum vacuum, the energy-density of the mediating vacuum field is Doppler-shifted by the relative cosmological recession. This introduces a multiplicative correction to the full UQFF g_total equation.

The Sombrero Galaxy (M104) is a nearby edge-on Sa-type spiral galaxy at a recession velocity of approximately +1890 km/s, corresponding to a spectroscopic heliocentric redshift of z = +0.0063. Its recession — the opposite of Andromeda's blueshift approach — provides the ideal test case for the positive-z branch of the universal κ(z) function.

This paper establishes:
1. The analytic form of κ_recession for z > 0.
2. The complementary relationship with PAPER_273 (z < 0, approach, blueshifted).
3. The universal UQFF Bidirectional Redshift Law valid for all real z > −1.
4. Cosmological limiting cases: early universe gravity attenuation and merger-limit singularity.

---

## 2. Theoretical Derivation

### 2.1 UQFF Vacuum Field Energy Under Recession

In UQFF, the gravitational term g_UQFF is proportional to the local vacuum energy density ρ_vac at the point of computation. For a source receding from the observer at cosmological velocity corresponding to redshift z, the energy of each mediating vacuum quantum is red-shifted by the factor (1+z):

$$E_{\text{obs}} = \frac{E_{\text{emit}}}{1+z}$$

Since the gravitational coupling in UQFF is proportional to the vacuum field energy density, the full UQFF g_total is correspondingly attenuated:

$$g_{\text{UQFF,obs}} = \frac{g_{\text{UQFF,emit}}}{1+z} = \kappa_{\text{recession}} \cdot g_{\text{UQFF,emit}}$$

Defining the recession damping factor:

$$\boxed{\kappa_{\text{recession}} = \frac{1}{1+z}}$$

For the Sombrero Galaxy:

$$\kappa_{\text{recession}}^{\text{Sombrero}} = \frac{1}{1 + 0.0063} = \frac{1}{1.0063} = 0.99374$$

### 2.2 Position in the computeG() Equation

The κ_recession factor enters the Sombrero UQFF Master Gravity Equation as an outer multiplier:

$$g_{\text{total}}(r, t) = \left[\, g_{\text{grav}} + U_{g,\text{sum}}^{(26)} + \Lambda_{\text{term}} + g_{\text{quantum}} + g_{\text{Lorentz}} + g_{\text{fluid}} + F_{\text{ring}}(t) + g_{\text{BH}} + g_{\text{exp}} + a_{\text{dust}} + g_{\text{DM}} \,\right] \cdot \kappa_{\text{recession}} \cdot \sigma_{\text{SC}}$$

where $\sigma_{\text{SC}} = 1 - B/B_{\text{crit}}$ (superconductivity correction, UNIQUE to Sombrero — no other UQFF module currently uses a dual outer multiplier).

**Gravitational attenuation for Sombrero:**

$$\Delta g_{\text{recession}} = g_{\text{UQFF}} \left(1 - \frac{1}{1.0063}\right) = 0.00626 \cdot g_{\text{UQFF}}$$

With g_base = G·M/r² = 2.382×10⁻¹⁰ m/s² and pre_sum_Ug ≈ 52·g_base ≈ 1.238×10⁻⁸ m/s²:

$$\Delta g \approx 0.00626 \times 1.238 \times 10^{-8} \approx 7.75 \times 10^{-11}\ \text{m/s}^2$$

---

## 3. Universal UQFF Bidirectional Redshift Law

Combining PAPER_273 (z = −0.001, Andromeda) and PAPER_277 (z = +0.0063, Sombrero), we arrive at the general statement:

**Universal UQFF Bidirectional Redshift Law:**

$$\kappa(z) = \frac{1}{1+z}, \quad z \in (-1, +\infty)$$

| Regime | z range | κ(z) | Physical meaning |
|--------|---------|------|-----------------|
| Approach/Blueshift | z < 0 | > 1 | UQFF gravity amplified (PAPER_273) |
| Rest frame | z = 0 | 1.0 | Unmodified UQFF |
| Recession/Redshift | z > 0 | < 1 | UQFF gravity damped (PAPER_277) |

### 3.1 κ(z) at Representative Redshifts

| Redshift z | κ(z) | Astrophysical Context |
|-----------|------|----------------------|
| z = −0.001 | 1.001001 | Andromeda (PAPER_273, approach) |
| z = 0 | 1.000000 | Rest frame / MW local group |
| z = +0.0063 | 0.99374 | **Sombrero M104 (PAPER_277)** |
| z = +0.1 | 0.90909 | Virgo Cluster periphery |
| z = +0.5 | 0.66667 | Intermediate cosmological |
| z = +1.0 | 0.50000 | Halfway epoch (t ≈ 5.8 Gyr) |
| z = +3.5 | 0.22222 | Epoch of reionisation |
| z → ∞ | → 0 | Big Bang limit — gravity switchoff |
| z → −1 | → ∞ | Singularity — merger/coalescence |

### 3.2 Cosmological Limiting Cases

**Early Universe (z → ∞):**

$$\lim_{z \to \infty} \kappa(z) = 0$$

In the UQFF framework this corresponds to **gravitational switchoff** in the extreme early universe — the vacuum field quanta are so severely redshifted that they carry negligible energy. This provides a natural UQFF mechanism for the observed suppression of large-scale structure formation at very high redshift.

**Merger Singularity (z → −1):**

$$\lim_{z \to -1} \kappa(z) = +\infty$$

For a source approaching the observer at the speed of light (z → −1 in the blueshift convention), the UQFF gravity diverges. This represents the merger or coalescence singularity, where the UQFF vacuum field collapses to a single point and gravitational focusing becomes unbounded.

---

## 4. Module Implementation

```cpp
// PAPER_277 — SOMBRERO_UQFF_MODULE.cpp, updateCache()
kappa_recession = 1.0 / (1.0 + z);   // z = +0.0063 → kappa_recession = 0.99374

// Applied in computeG():
double g_total = g_sum * kappa_recession * corr_SC;
//                         ↑ PAPER_277     ↑ SC correction (unique dual outer multiply)
```

**Computed values for Sombrero M104:**
- z = +0.0063
- kappa_recession = 0.99374
- Gravitational damping: −0.626% of g_UQFF (recession)

---

## 5. Key Constants Introduced

| Symbol | Value | Units | Description |
|--------|-------|-------|-------------|
| κ_recession | 0.99374 | dimensionless | UQFF recession damping factor for z=+0.0063 |
| z_Sombrero | +0.0063 | dimensionless | Heliocentric spectroscopic recession redshift |
| Δg_recession | ~7.75×10⁻¹¹ | m/s² | Absolute gravitational attenuation |

---

## 6. Physical Significance

1. **Bidirectional completeness:** PAPER_273 and PAPER_277 together show that the UQFF κ(z) correction is a universal single-parameter function. No new free parameters are introduced.

2. **Observable consequence for Sombrero:** The 0.626% attenuation is below current measurement precision for M104's rotation curves (~1–5% estimated observational errors), but is computable and constitutes a UQFF first-principles prediction.

3. **Cosmological implications:** The κ(z) law predicts that UQFF gravity was weaker in the early universe by the factor 1/(1+z), offering a potential explanation for the delayed gravitational collapse epoch without invoking modified gravity theories.

4. **Dual outer multiplier uniqueness:** Sombrero is the first UQFF module to employ *two* outer multipliers simultaneously (κ_recession × σ_SC), whose combined effect is:
   $$\kappa_{\text{recession}} \times \sigma_{\text{SC}} = 0.99374 \times (1 - 10^{-20}) \approx 0.99374$$

---

## 7. References

- PAPER_273: UQFF Gravitational Blueshift Amplifier for z < 0 (Andromeda Galaxy)
- PAPER_276: Friedmann H(z) Expansion Coupling in UQFF (Andromeda)
- SOMBRERO_UQFF_MODULE.cpp (UQFF 2.0, Session 77)
- Alpher, R. A., & Herman, R. C. (1948). *Physical Review*, 74, 1737. (Cosmological expansion)
- Ford, H. C. et al. (1996). ApJ, 458, 132. (Sombrero BH mass from gas kinematics)
- Emsellem, E. et al. (2004). MNRAS, 352, 721. (M104 stellar kinematics and redshift)

---

*UQFF 2.0 — All physics is additive. The κ_recession factor does not replace any prior term — it is a cosmological outer-multiplier consistent with vacuum field energy propagation theory. — Daniel T. Murphy, Session 77, March 2026.*
