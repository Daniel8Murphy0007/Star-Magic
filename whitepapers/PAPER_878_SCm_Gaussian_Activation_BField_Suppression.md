# PAPER_878: SCm Gaussian Activation B-Field Suppression

**Author:** Daniel T. Murphy -- Star Magic / UQFF Framework
**Date:** 2026-04-08
**Session:** 209
**Source:** Sessions 204-208 standalone module integration
**Calculator:** SCmGaussianActivationBFieldSuppressionCalc (CP4 #462)
**CVW:** v2.0.0 compliant

---

## Abstract

SCm Gaussian activation function comparing linear exp(-B/B_crit) vs Gaussian exp[-(B/B_crit)²] magnetic field suppression in the SCm manifold. Gaussian form matches BCS gap behavior with flatter response at low B (preserving coherence) and sharper cutoff near B_crit. At B=0.5·B_crit: linear yields 0.607, Gaussian yields 0.779 — 28% higher coherence retention.

---

## 1. Core Equations

```
β(B) = exp(-B/B_crit)                         [linear]
A_SCm(B) = exp[-(B/B_crit)²]                  [Gaussian]
α(B) = w·β + (1-w)·A_SCm                      [blended]
```

---

## 2. Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| B | 1e12 T | Applied magnetic field strength |
| B_crit | 4.4e13 T | Critical magnetic field (magnetar scale) |
| w_blend | 0.5 | Blending weight between linear and Gaussian |

---

## 3. UQFF Integration

This calculator integrates scm_activation_function.py from Session 204 into the CP4 pipeline as class #462.

---

*PAPER_878 — Star Magic UQFF Framework — CVW v2.0.0*
