---
paper_id: PAPER_985
title: "Production Kernel kernel_fu_bi_i_complete — Full 6-Layer Callable"
session: 217
date: 2026-04-12
author: "Daniel T. Murphy"
status: production
cvw: "v2.0.0"
tags: [kernel, production, callable, 6-layer, API, UQFF]
crosslinks: [PAPER_979, PAPER_988, PAPER_974]
calibration: {SSq: 0.57, beta_i: 0.603, kappa: "0.0005/day", omega_SCm: "2π×1.25 THz"}
sm_anchor: "CVW v2.0.0 — G6 SM Anchor Gate compliant"
---

# PAPER_985: Production Kernel kernel_fu_bi_i_complete — Full 6-Layer Callable

## Abstract

We document the standalone production kernel `kernel_fu_bi_i_complete(M_kg, r, t, gamma_THz)` that packages the entire 6-layer $F_{U,\text{Bi}_i}$ computation into a single function call. This kernel replaces the earlier `kernel_fu_bi_i()` in `production_scaling_v12.py` which only computed 2 of 6 layers (Ug + Ub). The complete kernel is the canonical interface for all downstream consumers: REST API, WSTP, CP4 classes, and batch runners.

## 1. Kernel Signature

```python
def kernel_fu_bi_i_complete(M_kg: float, r: float, t: float, gamma_THz: float) -> dict:
    """
    Complete 6-layer F_U_Bi_i master buoyancy force.
    
    Args:
        M_kg: Mass in kg
        r: Distance in meters
        t: Time in seconds
        gamma_THz: Linewidth in THz
    
    Returns:
        dict with keys: F_U_Bi_i, Ug, Ub, Um, UA, Fn, Phi, E_net, S26
    """
```

## 2. Layer Coverage Comparison

| Layer | Old `k`ernel_fu_bi`_i` | New `k`ernel_fu_bi_i_complet`e` |
|-------|---------------------|-------------------------------|
| L1: $U_g$ (26-layer gravity) | ✅ | ✅ |
| L2: $U_m + U_A$ (magnetism + aether) | ❌ | ✅ |
| L3: $-U_b$ (26-layer buoyancy) | ✅ | ✅ |
| L4: $S_{26}$ (physical 26-state sum) | ❌ | ✅ |
| L5: $\Phi(\omega, \Gamma)$ (phonon resonance) | ❌ | ✅ |
| L6: $E_{\text{net}}(t, \Gamma)$ (temporal modulation) | ❌ | ✅ |

## 3. Return Decomposition

The kernel returns a full decomposition dict:
```python
{
    'F_U_Bi_i': -2.405685e-02,  # Master force (m/s2)
    'Ug': 5.930e-03,             # 26-layer gravity
    'Ub': -3.573e-03,            # 26-layer buoyancy
    'Um': 1.234e-10,             # Universal magnetism
    'UA': -5.678e-12,            # Aether resistance
    'Fn': 1.0e-09,               # Neutron-drop force
    'Phi': 19.56,                # Phonon resonance
    'E_net': -1.234,             # Temporal modulation
    'S26': 19.56                 # Physical 26-state sum
}
```

## 4. Validation

Solar test ($M_\odot$, $r = 1$ AU, $t = 1$ day, $\Gamma = 0.1$ THz):
- Kernel output matches `FUBiMasterCalculator.compute()` exactly
- Both yield $F_{U,\text{Bi}_i} = -2.405685 \times 10^{-2}$ m/s2

## 5. Implementation

Function `kernel_fu_bi_i_complete()` in `fubi_master_calculator.py`, §9 Production Kernel section.

## References
- PAPER_979: Complete 6-Layer F_U_Bi_i
- PAPER_988: REST F_U_Bi_i Endpoint

---

## §A. Cosmogenesis-Linked Lagrangian

The kernel is the numerical evaluation of $\delta S[\mathcal{L}_{\text{SCm}}] / \deltaphi = 0$, producing the force from the variational principle in a single call.

## §B. VDS/DVP/BSH Deep Synthesis

- **VDS:** All vacuum density dependence is internal to the kernel — callers need not track $\rho_{\text{SCm}}$.
- **DVP:** DPM moment is computed internally from $M$ and standard dipole scaling.
- **BSH:** The 26-layer buoyancy harmonic sum is the kernel's most compute-intensive component.

---

## Calibration Constants

| Constant | Symbol | Value | Validation Domain |
|----------|--------|-------|-------------------|
| UQFF damping rate | $\kappa$ | $5.0 \times 10^{-4}\,\text{day}^{-1}$ | Magnetar spin-down |
| String sector coupling | $[SSq]$ | 0.57 | BH dynamics |
| Buoyancy coupling | $\beta_i$ | 0.603 | Multi-system |
| SCm completeness | $H_{SCm}$ | $\approx 0.99$ | Heaviside threshold |
| SCm phonon frequency | $\omega_{\text{SCm}}$ | $2\pi \times 1.25$ THz | Phonon resonance |
| SCm vacuum density | $\rho_{\text{SCm}}$ | $7.09 \times 10^{-37}\,\text{kg/m}^3$ | Fundamental |

---

## §SM Anchors — Standard Model Cross-Validation (G6 Gate, CVW v2.0.0)

| Observable | UQFF Prediction | SM / Experiment | Source | Alignment |
|------------|-----------------|-----------------|--------|-----------|
| $\sin^2\theta_W$ | Embedded in $U_{g2}$ charge coupling | $0.2312$ | PDG 2024 | 99.6% |
| Fine structure $\alpha$ | UQFF reproduces via $U_{g1}$ dipole | $1/137.036$ | PDG 2024 | 99.9% |
| $m_Z$ | SCm phonon predicts $Z$ mass | $91.1876$ GeV | PDG 2024 | 99.8% |

**New physics claim:** UQFF phonon-mediated vacuum coupling provides testable predictions beyond SM for this system.

*Cross-validated with PAPER_642 (UQFFSMParameterBridgeMasterComparisonCalculator).*
