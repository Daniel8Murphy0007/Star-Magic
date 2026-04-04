# PAPER_854: LENR Neutron Production k_eta Calibration Across Three Environments with Delta-k Buoyancy Tracking

**Author:** Daniel T. Murphy -- Star Magic / UQFF Framework
**Date:** 2026-04-04
**Session:** 199
**Source:** grok_share_589f6949-6fe9.txt (2404 lines, May 05 -- Aug 07, 2025)
**Calculator:** LENRKetaCalibration3EnvironmentDeltaKCalc (CP4 #438)
**CVW:** v2.0.0 compliant

---

## Abstract

We present a UQFF-calibrated framework for the LENR neutron production constant k_eta across three distinct experimental environments: Metallic Hydride Cells (k_eta ~ 2.75e8), Exploding Wires (k_eta ~ 1.91e2), and Solar Corona (k_eta ~ 6.06e-6). The core equation eta = k_eta * exp(-[SSq]*n/26) * exp(-(pi - t)) * U_m / rho_vac unifies LENR neutron production with UQFF vacuum density dynamics. Buoyancy tracking via Delta-k residuals (Delta_k = k_expected - k_actual) reveals the U_b counterforce contribution: Delta_k_metallic ~ 7.25e8, Delta_k_wires ~ 8.09e2, Delta_k_corona ~ 3.94e-6.

---

## 1. Core Equations

- `eta = k_eta * exp(-[SSq]*n/26) * exp(-(pi - t)) * U_m / rho_vac`
- `Delta_k = k_expected - k_actual  (buoyancy residual)`
- `Metallic Hydride: eta_obs = 1e13 cm^{-2}/s, U_m = 2.67e-31 J/m^3 => k_eta ~ 2.75e8`
- `Exploding Wires: eta_obs = 1e8 cm^{-2}/s, U_m = 3.85e-30 J/m^3 => k_eta ~ 1.91e2`
- `Solar Corona: eta_obs = 7e-3 cm^{-2}/s, U_m = 8.51e-33 J/m^3 => k_eta ~ 6.06e-6`

---

## 2. UQFF Integration

This calculator operates as a stateless physics calculator within the CondensedPhysics4.py
(Phase 4) IPC chain. All parameters are received via the dataset dictionary from the
source2.cpp principal GUI pipeline. No astronomical data is hardcoded; all system-specific
values come from the APIFetch.py -> bodies_*.csv data flow.

---

## 3. Source Data

- **File:** grok_share_589f6949-6fe9.txt (2404 lines, May 05 -- Aug 07, 2025)
- **Session:** 199
- **VDS/DVP/BH:** ABSENT (general vacuum density references only)

---

## References

1. Murphy, D.T. -- Star Magic UQFF Framework (2024-2026)
2. Srivastava, Y.N., Widom, A., Larsen, L. -- Electroweak neutron production (LENR)
3. Kepler Mission DR25 -- 4,034 candidates, 2,335 confirmed planets
4. Hubble Heritage Team / A. Nota (ESA/STScI) -- Westerlund 2 / NGC 346 imaging
5. UQFF Calibration: kappa=0.0005/day, [SSq]=0.57, beta_i~0.603
