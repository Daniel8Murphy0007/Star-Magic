# PAPER_899: Buoyancy Reversal Sign Flip Resonance

**Author:** Daniel T. Murphy -- Star Magic / UQFF Framework
**Date:** 2026-04-08
**Session:** 209
**Source:** Sessions 204-208 standalone module integration
**Calculator:** BuoyancyReversalSignFlipResonanceCalc (CP4 #483)
**CVW:** v2.0.0 compliant

---

## Abstract

Sweep of F_{U,Bi}/F_U from 0.1 to 0.9 detecting net_factor sign changes. The critical ratio R=0.5 marks the buoyancy-gravity balance point where net_factor = 0. Sign flips indicate expansion↔erosion phase transitions, analogous to cosmological phase changes.

---

## 1. Core Equations

```
net_factor(R) = 2R - 1
critical ratio at R = 0.5
sign flip = phase transition
```

---

## 2. Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| ratio_min | 0.1 | Minimum buoyancy ratio |
| ratio_max | 0.9 | Maximum buoyancy ratio |
| n_points | 9 | Number of sweep points |

---

## 3. UQFF Integration

This calculator integrates et_phonon_resonance.py from Session 208 into the CP4 pipeline as class #483.

---

*PAPER_899 — Star Magic UQFF Framework — CVW v2.0.0*
