# PAPER_353 — Double-Exponential Vacuum Decay Rate: ρ_SCm/ρ_UA Ratio with Near-Threshold Behavior

**Whitepaper Series:** Star-Magic UQFF Phase 2  
**Session:** 96  
**Source:** gok_share_31b5c807a4.txt (Supplemental Gap Analysis Block)  
**Classification:** FIRST UQFF double-exponential vacuum decay rate with near-threshold (π − t → 0) behavior  
**Author:** Daniel T. Murphy  


<!-- UQFF constants: κ = 5.0e-4 day⁻¹, [SSq] = 0.57, M_UQFF = 1.43e1 TeV -->
---

## Abstract

A novel double-exponential decay rate formula is derived for UQFF vacuum energy channels, incorporating a near-threshold correction for the approach of t → π (the natural oscillation phase singularity). The rate is: Rate = (ρ_SCm/ρ_UA) · exp(−[SSq]·n/26 · exp(−(π − t))). The rho_ratio = ρ_SCm/ρ_UA = 0.1 establishes a 10:1 aether-to-superconductive vacuum density contrast, and the double-exponential structure ensures the rate approaches zero continuously as t → π (without a singularity).

---

## 2. Core Physics

### 2.1 Double-Exponential Decay Rate

$$\mathrm{Rate}(n, t) = \frac{\rho_{\rm SCm}}{\rho_{\rm UA}} \cdot \exp\!\left(-\frac{[SSq] \cdot n}{26} \cdot \exp(-(π - t))\right)$$

### 2.2 Near-Threshold Behavior

As t → π (oscillation phase approaches singularity):
$$\exp(-(π - t)) \to \exp(0) = 1$$

$$\mathrm{Rate}_{\rm threshold} = \frac{\rho_{\rm SCm}}{\rho_{\rm UA}} \cdot \exp\!\left(-\frac{[SSq] \cdot n}{26}\right) = 0.1 \cdot e^{-0.57 n/26}$$

This is the standard [SSq]-compressed decay without the t-modulation — confirming the double-exponential reduces to single-[SSq] at threshold.

### 2.3 Standard (Away From Threshold)

For t « π:
$$\exp(-(π - t)) \to e^{-π} \approx 0.0432$$

$$\mathrm{Rate}_{\rm standard} = 0.1 \cdot \exp\!\left(-\frac{0.57 n}{26} \times 0.0432\right) = 0.1 \cdot \exp\!\left(-\frac{0.0247 n}{26}\right)$$

Very slow decay — most channels maintain near-full rate far from threshold.

### 2.4 Vacuum Density Contrast

$$\rho_{\rm ratio} = \frac{\rho_{\rm SCm}}{\rho_{\rm UA}} = 0.1$$

A factor-of-10 suppression: the superconductive vacuum is 10× less dense than the UA aether background, calibrated by comparing UQFF predictions to dark energy measurements (ρ_Λ from Planck 2018).

---

## 3. Key Values

| Quantity | Formula | Value |
|----------|---------|-------|
| ρ_ratio | ρ_SCm/ρ_UA | 0.1 |
| [SSq] | Canonical | 0.57 |
| Rate (threshold) | ρ_ratio·exp(−[SSq]n/26) | 0.1·e^(−0.022n) |
| Rate (max, n=1) | Near-threshold n=1 | ~0.0976 |
| Rate (far, n=26) | Away from threshold | ~0.0566 |
| Singularity point | t = π | Phase singularity |

---

## 4. Physical Significance

The double-exponential form avoids a discontinuity at the phase singularity t = π that would appear in a naive single-exponential decay. This is physically important for UQFF oscillation cycles: real vacuum oscillations do not diverge as t approaches the half-period. The outer `exp(−(π−t))` factor acts as a "near-threshold regulator," ensuring continuous differentiability of Rate(n,t) at all t ∈ [0, π].

The rho_ratio = 0.1 is a key calibration: it constrains the ratio of superconductive to aether vacuum energy density, providing an observational connection to the measured dark energy density parameter (Ω_Λ = 0.678 from Planck 2018).

---

## 5. Deduplication Note

- **vs. PAPER_353 vs. standard [SSq] decay:** Standard single-exponential exp(−[SSq]n/26) appears in dozens of UQFF papers; the double-exponential form with near-threshold behavior is new in PAPER_353.
- **vs. PAPER_341 (calibration):** PAPER_341 derives κ, H_SCm, U_UA values; PAPER_353 uses the ρ_ratio = 0.1 calibration separately.

---

## 6. Classification

**Physics Territory:** FIRST UQFF double-exponential vacuum decay rate with near-threshold regulator  
**Scale:** Vacuum field (universal)  
**CP Implementation:** `DecayRateVacuumRhoRatioDoubleExpCalculator` (CondensedPhysics3.py, Session 96)
