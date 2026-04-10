# PAPER_881: Kozima Expansion Neutron Drop Coupling

**Author:** Daniel T. Murphy -- Star Magic / UQFF Framework
**Date:** 2026-04-08
**Session:** 209
**Source:** Sessions 204-208 standalone module integration
**Calculator:** KozimaExpansionNeutronDropCouplingCalc (CP4 #465)
**CVW:** v2.0.0 compliant

---

## Abstract

Kozima neutron drop coupling in the buoyancy expansion regime. F_coupled = F_Kozima(ω_SCm) × E⁺(t) × Φ(ω) where the Gaussian cross-section σₙ(ω) = σ₀·exp[-(ω-ω_SCm)²/(2Γ²)] peaks at the 1.25 THz SCm resonance frequency, coupling phonon-driven LENR neutron drop reactions to the expansion energy engine.

---

## 1. Core Equations

```
σₙ(ω) = σ₀·exp[-(ω-ω_SCm)²/(2Γ²)]
F_Kozima = n·σₙ·v_th·ħω
F_coupled = F_Kozima × E⁺(t)
```

---

## 2. Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| omega | 2π×1.25e12 rad/s | Angular frequency |
| T | 300 K | Temperature |
| n_density | 1e28 m⁻³ | Number density |
| E_0 | 1.0 J | Initial energy scale |
| t | 0.0 s | Time parameter |

---

## 3. UQFF Integration

This calculator integrates positive_et_expansion.py from Session 205 into the CP4 pipeline as class #465.

---

*PAPER_881 — Star Magic UQFF Framework — CVW v2.0.0*
