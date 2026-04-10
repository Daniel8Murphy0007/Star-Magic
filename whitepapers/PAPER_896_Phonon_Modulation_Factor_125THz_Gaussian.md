# PAPER_896: Phonon Modulation Factor 1.25THz Gaussian

**Author:** Daniel T. Murphy -- Star Magic / UQFF Framework
**Date:** 2026-04-08
**Session:** 209
**Source:** Sessions 204-208 standalone module integration
**Calculator:** PhononModulationFactor125THzGaussianCalc (CP4 #480)
**CVW:** v2.0.0 compliant

---

## Abstract

Phonon modulation factor Φ_{1.25THz}(ω) = Φ₀·exp[-(ω-ω_SCm)²/(2Γ²)]·S₂₆([SSq]) with quality factor Q = ω_SCm/Γ = 6.25 (sharp resonance) and FWHM = 2Γ√(2ln2) ≈ 1.49 THz linewidth. Characterizes the SCm phonon resonance profile for LENR coupling.

---

## 1. Core Equations

```
Φ(ω) = Φ₀·exp[-(ω-ω_SCm)²/(2Γ²)]·S₂₆
Q = ω_SCm/Γ
FWHM = 2Γ√(2ln2)
```

---

## 2. Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| omega | 2π×1.25e12 rad/s | Angular frequency |

---

## 3. UQFF Integration

This calculator integrates et_phonon_resonance.py from Session 208 into the CP4 pipeline as class #480.

---

*PAPER_896 — Star Magic UQFF Framework — CVW v2.0.0*
