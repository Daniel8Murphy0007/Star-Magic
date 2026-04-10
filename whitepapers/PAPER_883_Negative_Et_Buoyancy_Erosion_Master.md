# PAPER_883: Negative E(t) Buoyancy Erosion Master

**Author:** Daniel T. Murphy -- Star Magic / UQFF Framework
**Date:** 2026-04-08
**Session:** 209
**Source:** Sessions 204-208 standalone module integration
**Calculator:** NegativeEtBuoyancyErosionMasterCalc (CP4 #467)
**CVW:** v2.0.0 compliant

---

## Abstract

Master equation for negative E(t) buoyancy-driven erosion. E⁻(t) = -E₀·exp(κt + [SSq]·t/26)·S₂₆·(1 - F_{U,Bi}/F_U). Symmetric counterpart to E⁺(t), active when gravity dominates buoyancy (ratio < 0.5). Applies to filament erosion, photoevaporation in M16, and GW damping scenarios.

---

## 1. Core Equations

```
E⁻(t) = -E₀·exp(κt + [SSq]·t/26)·S₂₆·(1 - R)
net_factor = 2R - 1
```

---

## 2. Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| E_0 | 1.0 J | Initial energy scale |
| t | 0.0 s | Time parameter |
| F_UBi_over_FU | 0.3 | Buoyancy-to-field ratio |

---

## 3. UQFF Integration

This calculator integrates negative_et_erosion.py from Session 205 into the CP4 pipeline as class #467.

---

*PAPER_883 — Star Magic UQFF Framework — CVW v2.0.0*
