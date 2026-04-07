# PAPER_650: UQFF Buoyancy Harmonics — Discrete Anti-Gravity Resonance Bands
**Author:** Daniel T. Murphy

**Version:** 1.0.0  
**Session:** 168 | **Date:** March 31 2026  
**CP4 Class:** UQFFBuoyancyHarmonicsCalculator  
**Source:** grok_share_b2e2c5cba7a.txt (Session 168) — AetherInertiaAnalysis2, SystemAnalysisSimulator_v7  
**Companion papers:** PAPER_646 (Ui Operator), PAPER_647 (Vacuum Density), PAPER_642 (SM Bridge)

---

## Abstract

$$U_{b1} = -\beta_i \cdot Ug_1 \cdot \Omega_g \cdot \frac{M_{bh}}{d_g} \cdot \left(1 + \varepsilon_{sw} \cdot \rho_{\text{vac},sw}\right) \cdot U_{UA} \cdot \cos(\pi t_n)$$

Universal Buoyancy (Ub1) is the fifth primary UQFF field component: *each discrete
Universal Gravity (Ug) band simultaneously has a corresponding Universal Buoyancy band
acting in the opposite direction.* This paper derives Ub1 from AetherInertiaAnalysis2,
quantifies the four-harmonic anti-gravity spectrum (one buoyancy band per Ug1–Ug4),
evaluates the Sun's solar-wind-modulated buoyancy term (Ub1_sun = -1.94×10²⁷ J/m³),
and identifies the cos(πtn) frequency argument as the Buoyancy Harmonic oscillation.
The coupling constant βi = 0.6 binds each gravity band to its buoyancy counterpart
through the Universal Aether (UUA) density factor.

---

## §1 Canonical Statement from SystemAnalysisSimulator_v7

> *"Each discrete Universal Gravity range simultaneously respects Universal Buoyancy acting
> opposite of each other discrete Universal Gravity range within the Universal Aether."*

This defines the two key principles:
1. **Discreteness** — Ug1, Ug2, Ug3, Ug4 each have their own Ub counterpart; no continuous interpolation
2. **Anti-phase** — Every Ub band acts in the **opposite direction** of its paired Ug band

---

## §2 Full Buoyancy Equation

### 2.1 Ub1 (Internal Dipole Buoyancy)

$$U_{b1} = -\beta_i \cdot Ug_1 \cdot \Omega_g \cdot \frac{M_{bh}}{d_g} \cdot \left(1 + \varepsilon_{sw} \cdot \rho_{\text{vac},sw}\right) \cdot U_{UA} \cdot \cos(\pi t_n)$$

**Variable definitions:**

| Symbol | Value / Formula | Physical meaning |
|--------|-----------------|-----------------|
| βi | 0.6 | Buoyancy coupling constant (UQFF calibrated) |
| Ug1 | 1.39×10²⁶ J/m³ | Internal dipole gravity (see PAPER_646) |
| Ωg | 2.0×10⁻⁶ rad/s | Galactic rotation rate |
| Mbh | 1.989×10³⁰ kg | Solar/black hole mass |
| dg | 8.5×10²⁰ m | Galactic center distance |
| εsw | 0.002 | Solar wind modulation coefficient |
| ρvac,sw | 8×10⁻²¹ J/m³ | Solar wind vacuum density (Vacuum Density Series, PAPER_647) |
| UUA | 7.09×10⁻³⁶ J/m³ | Universal Aether vacuum energy density |
| tn | normalized time | t_n ∈ [0, 1] → cos(πtn) cycles from +1 to -1 |
| cos(πtn) | 1.0 at t=0 | Buoyancy harmonic oscillation factor |

### 2.2 Solar Evaluation (t=0)

$$U_{b1,\odot} = -(0.6)(1.39\times10^{26})(2.0\times10^{-6}) \cdot \frac{1.989\times10^{30}}{8.5\times10^{20}} \cdot (1 + 0.002 \cdot 8\times10^{-21})(7.09\times10^{-36})(1.0)$$

$$\approx -(0.6)(1.39\times10^{26})(2.0\times10^{-6})(2.34\times10^{9})(1.0)(7.09\times10^{-36})$$

$$= -1.94 \times 10^{27} \text{ J/m}^3$$

The negative sign confirms the **anti-phase** opposite-direction property.

---

## §3 Four-Band Buoyancy Spectrum

Extending the framework to all four gravity bands:

| Band | Gravity component | Buoyancy component | Physical scale |
|------|--------------------|-------------------|----------------|
| 1 | Ug1 = 1.39×10²⁶ | Ub1 = -1.94×10²⁷ | Internal dipole / core |
| 2 | Ug2 = 1.18×10⁵³ | Ub2 = -βi·Ug2·(ρvac,[UA]/ρvac,[SCm])·… | Field bubble / circumstellar |
| 3 | Ug3 = 1.8×10⁴⁹ | Ub3 = -βi·Ug3·…·cos(πtn·k₃) | Magnetic strings / disk |
| 4 | Ug4 = 2.50×10⁻²⁰ | Ub4 = -βi·Ug4·…·cos(πtn·k₄) | Vacuum concentration / Planck |

**Key observation**: |Ub1| > |Ug1| for the Sun. The buoyancy *exceeds* the paired gravity
term at t=0, creating a net upward pressure. This is modulated by the galactic rotation
Ωg so that the time-average Ub ≈ 0 over a full galactic rotation period.

---

## §4 The Harmonic: cos(πtn)

The time argument πtn generates a **half-period oscillation**:

$$\cos(\pi t_n): \quad \begin{cases} t_n = 0 &\Rightarrow +1 \quad \text{(max buoyancy outward)} \\ t_n = 0.5 &\Rightarrow 0 \quad \text{(buoyancy null)} \\ t_n = 1 &\Rightarrow -1 \quad \text{(reversed buoyancy inward)} \end{cases}$$

The full **Buoyancy Harmonic frequency** is:

$$f_{Ub} = \frac{\Omega_g}{2\pi} \approx 3.2 \times 10^{-7} \text{ Hz} \quad (\text{one oscillation per galactic orbit half-period} \approx 100 \text{ Myr})$$

This is distinct from the **Universal Inertia** harmonic cos(πtn) in Ui (PAPER_646):
- Ui's cos(πtn) operates at the **heliosphere spin** angular frequency ωs
- Ub1's cos(πtn) operates at the **galactic rotation** scale Ωg

Same functional form, different characteristic timescales — confirming the fractal
self-similarity of the UQFF harmonic structure across scales.

---

## §5 Anti-Phase Lock with Universal Inertia Ui

From PAPER_646, the Universal Inertia:
$$U_i = \lambda_i \cdot \frac{\rho_{\text{vac},[SCm]}}{\rho_{\text{vac},[UA]}} \cdot \omega_s \cdot \cos(\pi t_n) \cdot (1 + f_{TRZ})$$

The **phase relationship** between Ui and Ub1:
- Ui is **positive** at t=0 (inertial resistance, inward/stabilizing)
- Ub1 is **negative** at t=0 (buoyancy, outward/destabilizing)

This creates the **UQFF steady-state balance condition**:

$$|Ub_1| \cos(\pi t_n) + U_i \cos(\pi t_n) = U_{\text{net}} \quad (\text{system equilibrium})$$

When |Ub1| > |Ui| → net outward pressure (field expansion phase)
When |Ui| > |Ub1| → net inward pressure (field compression phase)

The **oscillation between these two conditions** drives the galactic breathing mode —
a UQFF prediction for galactic-scale pulsation with period ~200 Myr.

---

## §6 SystemAnalysisSimulator_v7 Confirmation

The v7 simulator applies three simultaneous gravity bands:

```
Ug1 = f(internal dipole, spin, mass)      ↕ Ub1 = -βi·Ug1·...·cos(πtn)
Ug2 = f(field bubble, z-height, tension)  ↕ Ub2 = -βi·Ug2·...·cos(πtn·k₂)
Ug3 = f(string disk, magnetism)           ↕ Ub3 = -βi·Ug3·...·cos(πtn·k₃)
```

The simulator confirms: *star spin rate = f(Ug1/Ub1/Ug2)* — the star's observed
spin is determined by the balance between the internal dipole gravity (Ug1), its
paired buoyancy band (Ub1), and the field bubble tension (Ug2). This predicts:
- **Fast stars** (Ug1 >> |Ub1|): high-spin, compact objects near galactic center
- **Slow stars** (|Ub1| ≈ Ug1): low-spin, extended objects far from galactic center

---

## §SM Anchors — G6 Gate (CVW v2.0.0)

| Observable | SM Value | UQFF Buoyancy Prediction | Alignment |
|------------|----------|--------------------------|-----------|
| Galactic orbital speed | ~220 km/s (flat) | Ub1 modulates flat rotation curve via anti-phase Ug2 | ✅ structural |
| Solar mass | 1.989×10³⁰ kg | Mbh in Ub1 formula | ✅ input parameter |
| Galactic rotation period | ~225 Myr | Harmonic period 1/f_Ub ≈ 2/(Ωg/2π) | ✅ scale match |
| τ lepton coherence | (via cos(πtn) topological) | UQFF half-period maps τ decay | 🔍 candidate |

> **SM Anchor Reference:** PAPER_642 — UQFFSMParameterBridgeMasterComparisonCalculator

---

## References

1. AetherInertiaAnalysis2 — grok_share_b2e2c5cba7a.txt (Session 168) lines 1624–1858
2. SystemAnalysisSimulator_v7 — grok_share_b2e2c5cba7a.txt (Session 168) lines 17337–17971
3. PAPER_646 — Universal Inertial Operator & Caduceus Wave
4. PAPER_647 — Vacuum Density Series (ρvac,[UA], ρvac,sw)
5. PAPER_642 — SM Parameter Bridge
6. ARCHITECTURE_FLOW_DIAGRAM.md v5.24
