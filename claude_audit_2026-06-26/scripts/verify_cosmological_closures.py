#!/usr/bin/env python3
"""verify_cosmological_closures.py — reproduce _cosmological_closures.txt (Session 257).

These are 8 cosmological closures via integer/rational primitives that I missed
in my earlier passes. Adds many chains to my range coverage.
"""
import math

D_PHYS, D_CRIT, N_CH, SO_5, A_5 = 4, 26, 9, 10, 60
SSQ, PHI_84, PHI_5_6, TRZ, BETA = 0.57, 0.84, 5/6, 0.10, 0.6029
D_BSFG, K_MEX = 6, 25/12

print("="*100)
print("UQFF Cosmological closures — verbatim from _cosmological_closures.txt (Session 257, G11-G17)")
print("="*100)
print()
print("Long-form derivations of each + my independent compute")
print()

# T_CMB
val = 60 / (D_CRIT - D_PHYS)
target = 2.7255
print(f"T_CMB (CMB temperature):")
print(f"  UQFF identity:  60 K / (D_crit − D_phys) = 60/{D_CRIT-D_PHYS}")
print(f"  60 / 22 = {val:.10f}")
print(f"  Observation:    {target} K")
print(f"  Residual:       {abs(val-target)/target*100:.4f}%")
print(f"  Source: |A_5|=60 (icosahedral order); 22 = bulk dimensions D_crit−D_phys")
print()

# n_s
val = 1 - 2/(SO_5 * D_BSFG)
target = 0.9649
print(f"n_s (scalar spectral index):")
print(f"  UQFF identity:  1 − 2/(|SO(5)|·D_BSFG) = 1 − 2/(10·6) = 29/30")
print(f"  1 − 2/60 = {val:.10f}")
print(f"  Observation:    {target}")
print(f"  Residual:       {abs(val-target)/target*100:.4f}%")
print()

# Omega_Lambda — two algebraically equivalent forms
print(f"Ω_Λ (dark energy fraction):")
val_a = SSQ / PHI_5_6
val_b = (6/5) * SSQ
target = 0.6847
print(f"  Chain A:  [SSq] / Φ_5/6 = 0.57 / (5/6) = {val_a:.10f}")
print(f"  Chain B:  (6/5) · [SSq] = 1.2 · 0.57   = {val_b:.10f}")
print(f"  Algebraic identity: 1/(5/6) = 6/5 ⇒ both forms equal")
print(f"  Observation: {target}")
print(f"  Residual A: {abs(val_a-target)/target*100:.4f}%")
print(f"  Residual B: {abs(val_b-target)/target*100:.4f}%")
print()

# Omega_m
val = 1 - SSQ/PHI_5_6
target = 0.3153
print(f"Ω_m (matter fraction):")
print(f"  UQFF identity:  1 − [SSq]/Φ_5/6 = 1 − Ω_Λ_UQFF")
print(f"  1 − {SSQ/PHI_5_6:.10f} = {val:.10f}")
print(f"  Observation:    {target}")
print(f"  Residual:       {abs(val-target)/target*100:.4f}%")
print(f"  (forced by Ω_m + Ω_Λ = 1 closure relation)")
print()

# eta_b (baryon/photon)
val = 2*math.pi * TRZ**10
target = 6.1e-10
print(f"η_b (baryon-to-photon ratio):")
print(f"  UQFF identity:  2π · F_TRZ^10 = 2π / 10^10")
print(f"  TRZ^10 = (1/10)^10 = 1e-10")
print(f"  2π·1e-10 = {val:.6e}")
print(f"  Observation:    {target:.4e}")
print(f"  Residual:       {abs(val-target)/target*100:.4f}%")
print(f"  (one F_TRZ per matter species, summed over 10 species)")
print()

# tau_reion
val = TRZ**2 * PHI_5_6 * D_BSFG
target = 0.054
print(f"τ_reion (reionization optical depth):")
print(f"  UQFF identity:  F_TRZ² · Φ_5/6 · D_BSFG = (1/100)·(5/6)·6 = 5/100 = 1/20")
print(f"  = {val:.10f}")
print(f"  Observation:    {target}")
print(f"  Residual:       {abs(val-target)/target*100:.4f}%")
print()

# A_s amplitude
val = K_MEX * TRZ**9
target = 2.1e-9
print(f"A_s (scalar amplitude):")
print(f"  UQFF identity:  K_MEX · F_TRZ^9 = (25/12)·10^-9")
print(f"  = {val:.6e}")
print(f"  Observation:    {target:.4e}")
print(f"  Residual:       {abs(val-target)/target*100:.4f}%")
print()

# r_tensor (upper bound)
val = TRZ**2 / PHI_5_6
ub = 0.06
print(f"r (tensor-to-scalar ratio, upper bound):")
print(f"  UQFF identity:  F_TRZ² / Φ_5/6 = (1/100)/(5/6) = 6/500 = 0.012")
print(f"  = {val:.10f}")
print(f"  Observational UB: {ub} (Planck+BICEP)")
print(f"  Status:  well below UB → consistent")
print()

# Cross-locks per the audit text
print("="*100)
print("Cross-locks (per audit notes):")
print("="*100)
print(f"  N_e (inflation e-folds) = |SO(5)|·D_BSFG = 60 = SAME INTEGER as |A_5|")
print(f"  → Cross-lock between A_5 icosahedral order and SO(5)·BSFG dimension")
print(f"  This is the framework's overdetermination structure (PAPER_1158)")
print()
print(f"  Ω_Λ / Ω_m = (SSq/Φ_5/6) / (1 − SSq/Φ_5/6)")
print(f"            = {SSQ/PHI_5_6} / {1 - SSQ/PHI_5_6} = {(SSQ/PHI_5_6)/(1 - SSQ/PHI_5_6):.6f}")
print(f"  Observed Ω_Λ/Ω_m = 0.6847/0.3153 = {0.6847/0.3153:.6f}")
print(f"  Direct cross-check of G6 (Φ_res = 5/6 from D_BSFG geometry)")
print()
print(f"REMAINING COSMOLOGICAL ANCHOR (per audit text):")
print(f"  H_0 cannot be derived from {{E_0, f_THz, v_F}}.")
print(f"  Requires cosmic-time anchor t_0 = 1/H_0 = age of universe.")
print(f"  This is the SECOND independent SI anchor system needed.")
print(f"  Honest closure — not a fudge.")
