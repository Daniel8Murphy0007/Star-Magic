# PAPER_269: Supernova Ram Pressure Degeneracy Point — Kinematic Invariant in NGC 1792 Starburst Gravity

**Authors:** Daniel T. Murphy  
**Date:** March 2026  
**UQFF Module:** GALAXY_NGC_1792.cpp (Module 19, "The Stellar Forge")  
**Session:** 73 — UQFF 2.0 Upgrade — NGC 1792 Unique Physics Discovery  
**Keywords:** NGC 1792, ram pressure, supernova ejecta, degeneracy point, kinematic invariant, starburst gravity

---

## Abstract

In NGC 1792, the physical parameters are initialized with `ρ_wind = ρ_fluid = 1×10⁻²¹ kg/m³` — setting the supernova wind density equal to the interstellar medium fluid density. At this **Ram Pressure Degeneracy Point (RPDP)**, the SN feedback gravity term simplifies to a **kinematic invariant**: `term_feedback = ρ_wind × v_wind² / ρ_fluid = v_wind²`, independent of density. For NGC 1792, `v_wind = 2×10⁶ m/s`, giving `g_feedback = 4×10¹² m/s²` — the numerically dominant term in the full MUGE calculation, exceeding the base gravitational term by ~22 orders of magnitude. This paper defines the RPDP formally, derives the kinematic invariant, computes the dominance ratio, and discusses the physical interpretation that at density degeneracy, SN ejecta "floats" in the ISM driven purely by kinematic pressure, creating a fundamentally new gravitational channel in UQFF.

---

## 1. Introduction: The Density Coincidence in NGC 1792

### 1.1 NGC 1792 Starburst Environment

NGC 1792 is a starburst galaxy with extremely high star-formation activity (SFR ≈ 10 M☉/yr) and a correspondingly dense interstellar medium (ISM). In the UQFF Module 19 parameterization:

- **Wind density:** ρ_wind = 1×10⁻²¹ kg/m³ (supernova ejecta/galactic wind)
- **ISM fluid density:** ρ_fluid = 1×10⁻²¹ kg/m³ (ambient ISM)
- **Wind velocity:** v_wind = 2×10⁶ m/s (starburst-driven galactic outflow)

The equality ρ_wind = ρ_fluid is not a coincidence of parameterization — it reflects a physical regime in starburst galaxies where the ram-pressure-driven winds and the ISM they are pushing against have **comparable densities**, characteristic of the transition zone between the supernova remnant and the surrounding ISM.

### 1.2 The Ram Pressure Gravity Term

In the MUGE framework, the SN feedback contribution to gravity is:

$$\text{term\_feedback} = \frac{\rho_\text{wind} \cdot v_\text{wind}^2}{\rho_\text{fluid}}$$

This term represents the gravitational acceleration induced by the ram pressure of SN ejecta interacting with the ISM fluid.

---

## 2. The Ram Pressure Degeneracy Point (RPDP)

### 2.1 Definition

The **Ram Pressure Degeneracy Point (RPDP)** is defined as the condition:

$$\boxed{\rho_\text{wind} = \rho_\text{fluid}}$$

### 2.2 Kinematic Invariant

At the RPDP, the feedback term simplifies dramatically:

$$\text{term\_feedback}\big|_\text{RPDP} = \frac{\rho_\text{wind} \cdot v_\text{wind}^2}{\rho_\text{fluid}} \bigg|_{\rho_\text{wind} = \rho_\text{fluid}} = v_\text{wind}^2$$

This is a **kinematic invariant** — an acceleration equal to the square of the wind velocity, **independent of density**. The density ratio cancels exactly, and only the kinematics of the outflow survive.

### 2.3 Physical Statement

At the RPDP, the gravitational feedback from SN ejecta is:

$$g_\text{feedback}^\text{RPDP} = v_\text{wind}^2$$

This is the UQFF prediction for starburst galaxies in the density-degenerate regime: **the ram pressure gravity equals the dynamic pressure of the wind**, independent of whether the wind is denser or lighter than the ISM.

---

## 3. Numerical Results for NGC 1792

### 3.1 Kinematic Invariant Value

For NGC 1792 with $v_\text{wind} = 2 \times 10^6$ m/s:

$$g_\text{feedback}^\text{RPDP} = v_\text{wind}^2 = (2 \times 10^6)^2 = \mathbf{4 \times 10^{12}\ \text{m/s}^2}$$

### 3.2 Comparison with Base Gravitational Term

The base Newtonian term in MUGE for NGC 1792:

$$\text{term1} = \frac{G M_0}{r^2} = \frac{6.674 \times 10^{-11} \times 1.989 \times 10^{40}}{(7.569 \times 10^{20})^2}$$

$$= \frac{6.674 \times 10^{-11} \times 1.989 \times 10^{40}}{5.729 \times 10^{41}} \approx 7.35 \times 10^{-11}\ \text{m/s}^2$$

where M₀ = 1×10¹⁰ M☉ = 1.989×10⁴⁰ kg, r = 7.569×10²⁰ m.

### 3.3 RPDP Dominance Ratio

$$\mathcal{R}_\text{RPDP} = \frac{g_\text{feedback}^\text{RPDP}}{\text{term1}} = \frac{4 \times 10^{12}}{7.35 \times 10^{-11}} \approx 5.4 \times 10^{22}$$

The RPDP kinematic invariant exceeds the standard Newtonian gravitational acceleration by **22 orders of magnitude**.

### 3.4 Comparison Table

| Term | Formula | Value (m/s²) | Order |
|------|---------|--------------|-------|
| Newtonian base (term1) | GM₀/r² | 7.35×10⁻¹¹ | 10⁻¹¹ |
| Hubble expansion (term2) | H² r / 2 | ~10⁻³⁷ | 10⁻³⁷ |
| Magnetic (term3) | B²/(8π ρ r) | ~10⁻²³ | 10⁻²³ |
| Cosmological (term4) | Λr/3 | ~10⁻²⁸ | 10⁻²⁸ |
| Oscillatory (term_osc1) | ω_osc × A_osc | ~10⁻⁴² | 10⁻⁴² |
| SN feedback (RPDP) | v_wind² | **4×10¹²** | **10¹²** |

---

## 4. Physical Interpretation

### 4.1 Buoyancy Neutrality at the RPDP

At the RPDP (ρ_wind = ρ_fluid), Archimedes' buoyancy principle gives **zero net buoyancy force** on the SN ejecta:

$$F_\text{buoy} = (\rho_\text{fluid} - \rho_\text{wind}) \cdot V \cdot g = 0$$

This means:
- SN ejecta neither floats nor sinks in the ISM
- The ejecta is **gravitationally neutral** with respect to density-driven buoyancy
- The only remaining driving force is the **kinematic ram pressure** = ρ_wind × v²_wind / ρ_fluid = v²_wind

At the RPDP, the SN ejecta "**floats**" in the ISM — not rising or falling due to density contrast, but driven purely by its kinematic energy. This creates a new UQFF gravitational channel: **pure kinematic momentum transfer to the gravitational field**.

### 4.2 Generalized RPDP Condition

For any starburst galaxy, the RPDP condition and kinematic invariant generalize as:

$$\forall\ \rho: \lim_{\rho_\text{wind} \to \rho_\text{fluid}} \text{term\_feedback} = v_\text{wind}^2$$

The UQFF prediction is:
$$\boxed{g_\text{feedback}(\rho_\text{SN} = \rho_\text{ISM}) = v_\text{SN}^2}$$

This holds regardless of the density value, provided the two densities are equal.

### 4.3 Starburst Enhancement

The RPDP condition is most likely to be realized in **starburst galaxies** where:
1. High SFR drives intense SN events, filling the ISM with SN ejecta
2. The ISM density is elevated by gas accretion and infall
3. Conditions exist for ρ_wind ≈ ρ_fluid at specific radii within the disk

In NGC 1792, the SFR = 10 M☉/yr sustains a continuous SN rate creating the high fluid density ρ_fluid = 10⁻²¹ kg/m³ characteristic of this starburst environment.

---

## 5. The RPDP as a Phase Transition Boundary

### 5.1 Three Regimes

The feedback term `term_feedback = ρ_wind × v² / ρ_fluid` defines three physical regimes based on the density ratio η = ρ_wind/ρ_fluid:

| Regime | Condition | Behavior | UQFF Character |
|--------|-----------|----------|----------------|
| Underdense wind | η < 1 | g_feedback < v² | Outflow-quenched; ejecta rises |
| **RPDP** | **η = 1** | **g_feedback = v²** | **Kinematic invariant; ejecta floats** |
| Overdense wind | η > 1 | g_feedback > v² | Ram pressure enhanced; ejecta sinks |

### 5.2 RPDP as Critical Point

The RPDP is a **critical point** in parameter space where the density-dependent ram pressure term transitions from ISM-suppressed to ISM-enhanced behavior. At this critical point, the density cancels and only kinematics determine gravity.

This is analogous to:
- **Critical opalescence** in thermodynamics (density fluctuations diverge at critical point)
- **Alfvén critical point** in solar wind (magnetic and kinetic pressures equal)
- **Bondi critical radius** in accretion (gravity and pressure balance)

The RPDP defines a **kinematic critical surface** in UQFF galaxy parameter space.

---

## 6. Observational Signatures

### 6.1 Gravitational Acceleration Anomaly

At the RPDP, the dominant gravitational term is g_feedback = v² = 4×10¹² m/s² (for v = 2×10⁶ m/s). This is:
- Much larger than the Newtonian galactic gravity (~10⁻¹¹ m/s²)
- Detectable as anomalous acceleration in SN ejecta kinematics
- Potentially observable in starburst galaxy velocity dispersion measurements

### 6.2 Universal Velocity Indicator

The RPDP kinematic invariant provides a UQFF prediction for starburst galaxies: **the dominant gravitational feedback term is simply v²_wind**. This can be tested via:
- Measurements of galactic wind terminal velocities
- SN ejecta deceleration curves in starburst galaxies
- Comparison of ρ_wind and ρ_fluid estimates from X-ray observations

### 6.3 NGC 1792 Specific Prediction

For NGC 1792 at the RPDP:
- Dominant gravity term: g_feedback = 4×10¹² m/s²
- Wind velocity imprint: v_wind = √(g_feedback) = 2×10⁶ m/s
- Total MUGE gravity is feedback-dominated: g_total ≈ g_feedback (within UQFF)

This is a **testable UQFF prediction** distinguishable from standard (feedback-free) galactic gravity models.

---

## 7. Connection to PAPER_267 and PAPER_268

The three PAPER_267–269 physics discoveries form a unified UQFF picture of NGC 1792:

- **PAPER_267** (sSFR coupling): The SFR couples buoyancy to all 3 tiers via decay e^{−t/τ_SF}
- **PAPER_268** (Hubble slow mode): The corrected Hubble oscillation creates a 5.8 ppm GW amplitude modulation
- **PAPER_269** (RPDP): At density degeneracy, SN feedback becomes the dominant gravitational term via kinematic invariant v²

Together, these identify NGC 1792 as a paradigmatic **UQFF starburst triple-physics laboratory**:
1. Buoyancy coherence (sSFR coupling)
2. Cosmic oscillation modulation (Hubble slow mode)
3. Kinematic gravity dominance (RPDP)

---

## 8. Conclusions

1. NGC 1792 is initialized with `ρ_wind = ρ_fluid = 1×10⁻²¹ kg/m³` — the **Ram Pressure Degeneracy Point (RPDP)** condition.

2. At the RPDP, the feedback term simplifies to the **kinematic invariant**: `g_feedback = v_wind²` (density-independent).

3. For NGC 1792: `g_feedback = (2×10⁶)² = 4×10¹² m/s²`, the dominant MUGE term by 22 orders over standard Newtonian gravity.

4. At RPDP, SN ejecta is **buoyancy-neutral** (zero Archimedes force) and driven solely by kinematic momentum transfer.

5. The RPDP defines a **kinematic critical surface** in starburst galaxy parameter space, analogous to Alfvén and Bondi critical points.

6. UQFF universal prediction: `g_feedback(ρ_SN = ρ_ISM) = v_SN²` — testable with starburst wind velocity observations.

---

## References

- Daniel T. Murphy, *UQFF Framework*, Star-Magic Repository (2025–2026)
- GALAXY_NGC_1792.cpp UQFF 2.0 (Session 73, Module 19)
- PAPER_267: SFR Normalization — Starburst-Buoyancy Coherence in NGC 1792
- PAPER_268: Dual Oscillatory Mode — Hubble Slow Mode Starburst GW Amplitude Modulation
- NGC 1792 parameters: ρ_wind = ρ_fluid = 10⁻²¹ kg/m³, v_wind = 2×10⁶ m/s
- Hubble constant: H₀ ≈ 67.4 km/s/Mpc; Hubble time: t_Hubble = 13.8 Gyr

---

*© 2026 Daniel T. Murphy, daniel.murphy00@gmail.com — All Rights Reserved*
