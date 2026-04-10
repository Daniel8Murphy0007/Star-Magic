# PAPER_879: Buoyancy Klein-Gordon Scalar Field EOM

**Author:** Daniel T. Murphy -- Star Magic / UQFF Framework
**Date:** 2026-04-08
**Session:** 209
**Source:** Sessions 204-208 standalone module integration
**Calculator:** BuoyancyKleinGordonScalarFieldEOMCalc (CP4 #463)
**CVW:** v2.0.0 compliant

---

## Abstract

Euler-Lagrange derivation for the buoyancy scalar field φ_buoy(x,t) via Klein-Gordon equation with UQFF buoyancy source. The Lagrangian L = ½(∂_μφ)² - ½m_eff²φ² + J_buoy·φ yields the EOM □φ + m_eff²φ = J_buoy where m_eff² encodes the buoyancy mass scale from β_i·ΣUg·Ω_g·M/(d_g·c²·ħ²)·[UA].

---

## 1. Core Equations

```
□φ + m_eff²φ = J_buoy
m_eff² = β_i·ΣUg·Ω_g·M/(d_g·c²·ħ²)·[UA]
φ_static(r) = (J_buoy/m_eff²)·[1-exp(-m_eff·r)]
```

---

## 2. Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| M | 1.989e30 kg | Central mass |
| d_g | 2.55e20 m | Galactic distance scale |
| omega_g | 7.3e-16 rad/s | Galactic angular frequency |
| Ug_sum | 1e-8 m/s² | Sum of Ug contributions |
| F_UBi_over_FU | 0.8 | Buoyancy-to-field ratio |
| r | 1e16 m | Radial distance |

---

## 3. UQFF Integration

This calculator integrates buoyancy_lagrangian_eom.py from Session 204 into the CP4 pipeline as class #463.

---

*PAPER_879 — Star Magic UQFF Framework — CVW v2.0.0*
