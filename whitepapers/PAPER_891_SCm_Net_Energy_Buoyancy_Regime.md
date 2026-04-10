# PAPER_891: SCm Net Energy Buoyancy Regime

**Author:** Daniel T. Murphy -- Star Magic / UQFF Framework
**Date:** 2026-04-08
**Session:** 209
**Source:** Sessions 204-208 standalone module integration
**Calculator:** SCmNetEnergyBuoyancyRegimeCalc (CP4 #475)
**CVW:** v2.0.0 compliant

---

## Abstract

Net energy in the SCm vacuum buoyancy regime. E_net,SCm = ρ_SCm(t)·V·c²·(2R-1). Determines whether the system is in an expansion regime (nebulae, cosmogenesis) or erosion regime (filaments, cavities) based on the buoyancy ratio R = F_{U,Bi}/F_U.

---

## 1. Core Equations

```
E_net,SCm = ρ_SCm(t)·V·c²·(2R-1)
```

---

## 2. Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| t | 0.0 s | Time parameter |
| V | 1e48 m³ | Volume |
| F_UBi_over_FU | 0.8 | Buoyancy-to-field ratio |

---

## 3. UQFF Integration

This calculator integrates et_scm_vacuum.py from Session 207 into the CP4 pipeline as class #475.

---

*PAPER_891 — Star Magic UQFF Framework — CVW v2.0.0*
