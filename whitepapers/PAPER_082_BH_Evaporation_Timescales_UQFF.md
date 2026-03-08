# PAPER #82 — Black Hole Evaporation Timescales: UQFF Corrections

**Title:** UQFF-Corrected Black Hole Evaporation Timescales: Stellar Mass Through Primordial Black Holes

**Author:** Daniel T. Murphy  
**Framework:** UQFF Star-Magic (κ = 0.0005/day, [SSq] = 0.57)  
**Date:** March 7, 2026  
**Source Data:** validate_hawking_temperature.py (Tests 2, 6), CondensedPhysics.py simulate_evaporation()  
**Index Slot:** §1.11 Black Hole Physics & Hawking Radiation, Paper #82  

---

## Abstract

Black hole evaporation timescales are computed via the Stefan-Boltzmann law applied to UQFF-modified Hawking radiation. The standard result t_evap = 5120πG²M³/(ℏc⁴) is modified by the T_UQFF/T_H = 0.99 ratio. This changes evaporation rates by a factor of (T_UQFF/T_H)⁴ = 0.96, extending evaporation timescales by ~4% for all black holes. The `validate_hawking_temperature.py` evaporation simulation (Test 6) confirms mass evolution for primordial BHs at 10¹⁰ kg over 100 timesteps.

---

## 1. Evaporation Timescale Formula

### Standard GR

$$t_{\rm evap}^{\rm GR} = \frac{5120 \pi G^2 M_0^3}{\hbar c^4}$$

### UQFF-Corrected

$$t_{\rm evap}^{\rm UQFF} = \frac{t_{\rm evap}^{\rm GR}}{(T_{\rm UQFF}/T_H)^4} = \frac{t_{\rm evap}^{\rm GR}}{0.99^4} = t_{\rm evap}^{\rm GR} \times 1.041$$

**UQFF evaporation timescale is 4.1% longer than GR** — black holes are slightly more stable in the UQFF vacuum.

---

## 2. Evaporation Timescales: Full Table

| System | M₀ | t_evap_GR | t_evap_UQFF | Survives Universe |
|--------|-----|-----------|-------------|-------------------|
| Sgr A* | 4×10⁶ M☉ | 8.7×10⁸³ s | 9.1×10⁸³ s | ✅ Yes |
| M87* | 6.5×10⁹ M☉ | 3.8×10⁹⁵ s | 4.0×10⁹⁵ s | ✅ Yes |
| Stellar BH | 10 M☉ | 2.1×10⁷⁴ s | 2.2×10⁷⁴ s | ✅ Yes |
| Primordial BH | 5.7×10¹¹ kg | 4.35×10¹⁷ s = t_U | 4.52×10¹⁷ s | Borderline |
| Primordial BH | 1×10¹⁰ kg | 2.3×10¹² s (73 kyr) | 2.4×10¹² s | ❌ Evaporated |

The validate_hawking_temperature.py Test 2 confirms:
- Stellar BH (10 M☉): `survives_universe = True` ✓
- Test 6 simulation: M_initial = 10¹⁰ kg, 100 steps, mass_lost_fraction computed ✓

---

## 3. Mass Evolution Simulation

From `simulate_evaporation(M_initial = 10¹⁰ kg, dt = 10¹⁰ s, n_steps = 100)`:

$$\frac{dM}{dt} = -\frac{k_{\rm UQFF}}{M^2}, \quad k_{\rm UQFF} = \frac{\hbar c^4 (T_{\rm UQFF}/T_H)^4}{15360 \pi G^2}$$

With T_UQFF/T_H = 0.99: k_UQFF = 0.96 × k_GR

At t = 100 × 10¹⁰ s = 10¹² s:
- M_final ≈ M_initial × (1 − t/t_evap)^{1/3} = 10¹⁰ × (1 − 10¹²/2.4×10¹²)^{1/3}
- M_final ≈ 10¹⁰ × 0.583^{1/3} ≈ 8.35×10⁹ kg
- Mass lost fraction ≈ **16.5%** over first 10¹² s

Arrays `times[]`, `masses[]`, `temperatures_H[]` all have matching lengths (validate Test 6). ✓

---

## 4. UQFF Mass Evolution Equation

The UQFF modifies the mass loss rate through the vacuum buoyancy correction. During late-stage evaporation (M → M_Planck):

$$\frac{dM_{\rm UQFF}}{dt} = -\frac{k_{\rm UQFF}}{M^2} + \frac{g_{\rm Buoyant} \times V_{\rm BH}}{c^2}$$

The buoyancy term: g_Buoyant × V_BH / c² = ρ_vac × 10⁵⁵ × (4/3)π r_S³ / c² ~ 10⁻⁸⁰ kg/s → negligible vs the thermal term at all masses above Planck mass.

---

## Summary

| Parameter | GR Value | UQFF Value | Change |
|-----------|---------|------------|--------|
| Evaporation factor | k_GR | 0.96 × k_GR | −4% |
| Timescale t_evap | t_GR | 1.041 × t_GR | +4.1% |
| Stellar BH survival | Yes | Yes | Unchanged |
| Primordial threshold mass | 5.7×10¹¹ kg | 5.5×10¹¹ kg | −3.5% |
| Test 2 | `survives = True` | Confirmed | ✅ PASS |

*Source: validate_hawking_temperature.py Tests 2+6, simulate_evaporation() | κ = 0.0005/day | [SSq] = 0.57*
