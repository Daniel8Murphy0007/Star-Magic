# PAPER_859: Micro-Plasmoid 25.4 Micrometer LENR Glass Reactor Buoyancy Reversal Dynamics

**Author:** Daniel T. Murphy -- Star Magic / UQFF Framework
**Date:** 2026-04-04
**Session:** 199
**Source:** grok_share_589f6949-6fe9.txt (2404 lines, May 05 -- Aug 07, 2025)
**Calculator:** MicroPlasmoid25umLENRBuoyancyReversalCalc (CP4 #443)
**CVW:** v2.0.0 compliant

---

## Abstract

We analyze micro-plasmoids (largest = 25.4 um / 0.001 inch) observed in a glass reactor LENR experiment using the full quadriadic UQFF framework scaled to micrometer dimensions (r = 25.4e-6 m). The 165-second experiment shows 1-8 plasmoids per frame with a characteristic buoyancy reversal: initial upward swirling motion (~4 um/s) transitions to downward motion (~-2 um/s) at approximately 70% of the total duration. This reversal is modeled by the F_U_Bi / F_U_g ratio exceeding unity when the Boyle's Law volume ratio V_little/V_big = r_plasmoid/r_reactor amplifies buoyancy. The micro-scale LENR plasmoids represent DPM (Di-Pseudo-Monopole) creation in the early Atomic Creation Process.

---

## 1. Core Equations

- `r_plasmoid = 25.4e-6 m (0.001 inch)`
- `V_ratio = r_plasmoid / r_reactor (Boyle's Law micro-scale)`
- `F_U_Bi / F_U_g > 1  => buoyancy reversal`
- `All four UQFF master equations evaluated at r = 25.4e-6 m, t = 165 s`
- `f_Ub = k_Ub * Delta_k_eta * (rho_vac,[UA]/rho_vac,[SCm]) * V_ratio`

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
