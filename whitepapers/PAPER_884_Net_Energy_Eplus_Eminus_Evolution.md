# PAPER_884: Net Energy E⁺/E⁻ Evolution

**Author:** Daniel T. Murphy -- Star Magic / UQFF Framework
**Date:** 2026-04-08
**Session:** 209
**Source:** Sessions 204-208 standalone module integration
**Calculator:** NetEnergyEplusEminusEvolutionCalc (CP4 #468)
**CVW:** v2.0.0 compliant

---

## Abstract

Net energy evolution E_net(t) = E⁺(t) + E⁻(t) = E₀·exp(κt + [SSq]·t/26)·S₂₆·(2R-1). The sign is determined by the buoyancy ratio R: R>0.5 yields expanding universe, R<0.5 yields eroding/contracting regime, R=0.5 is perfectly balanced. Verifies the identity E⁺ + E⁻ ≡ E_net.

---

## 1. Core Equations

```
E_net = E⁺ + E⁻ = E₀·exp(...)·S₂₆·(2R-1)
verified identity
```

---

## 2. Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| E_0 | 1.0 J | Initial energy scale |
| t | 0.0 s | Time parameter |
| F_UBi_over_FU | 0.8 | Buoyancy-to-field ratio |

---

## 3. UQFF Integration

This calculator integrates negative_et_erosion.py from Session 205 into the CP4 pipeline as class #468.

---

*PAPER_884 — Star Magic UQFF Framework — CVW v2.0.0*
