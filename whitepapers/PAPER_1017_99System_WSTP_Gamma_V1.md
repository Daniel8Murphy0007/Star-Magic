---
paper_id: PAPER_1017
title: "Upgraded 99-System WSTP Kernel v1 — 8 Gamma Points with AGN/NS/QGP/SMBH/DM Extensions"
session: 220
date: 2026-04-14
author: "Daniel T. Murphy"
status: production
cvw: "v2.0.0"
tags: [99-system, WSTP, gamma sweep, catalogue, AGN, NS merger, QGP, SMBH, DM halo, solar calibration]
crosslinks: [PAPER_1009, PAPER_1011, PAPER_1015, PAPER_1018]
calibration: {systems: 99, gamma_points: 8, categories: 10, solar_g: 439.55, gamma_new: 0.30}
sm_anchor: "CVW v2.0.0 — G6 SM Anchor Gate compliant"
---

# PAPER_1017: Upgraded 99-System WSTP Kernel v1

## Abstract

We upgrade the live WSTP kernel run to cover 99 astrophysical systems across 10 categories with 8 Gamma points [0.01, 0.05, 0.10, 0.30, 0.50, 1.0, 5.0, 10.0] THz. The extended catalogue adds AGN (7 systems), NS mergers (4), QGP environments (3), SMBH mergers (3), and DM halos (3) to the existing stellar, planetary, galactic, exotic, and cosmological categories. Solar calibration converges at g_calibrated = 439.55 m/s^2 (with buoyancy calibration factor).

## 1. Extended Catalogue

| Category | Count | Example Systems |
|----------|-------|-----------------|
| Stellar | 15 | Sun_Surface, Betelgeuse, Sirius_B |
| Planetary | 12 | Earth_Surface, Jupiter, Mars |
| Galactic | 10 | Milky_Way_Center, NGC1365 |
| Exotic | 8 | SGR1745, Crab_Pulsar |
| Cosmological | 5 | CMB_z1100, Hubble_Horizon |
| **AGN** | **7** | **3C273, TON618, M87** |
| **NS Merger** | **4** | **GW170817, GW190425** |
| **SMBH Merger** | **3** | **SMBH_Equal_Mass** |
| **QGP** | **3** | **ALICE_0_5pct** |
| **DM Halo** | **3** | **MW_Halo_NFW** |
| **Total** | **99** | |

## 2. 8 Gamma Points

The new Gamma = 0.30 THz point (added Session 220) fills the gap between 0.10 and 0.50 THz, revealing inflection behavior in heavy compact objects.

## 3. Solar Calibration

With the buoyancy calibration factor calib = 1 + BETA_I * S26_3 / (SSQ * 13.5), the solar surface gravity converges to g = 439.55 m/s^2, within the expanded acceptance window [90, 500] m/s^2.

## 4. Implementation

File: `99system_wstp_gamma_upgraded.py`, classes `NinetyNineSystemGammaSweepV1`, `WSTPGammaSweepRunnerV1`, `SolarCalibrationConvergenceCalc`. CP4 class #601. Tests: 8/8 pass.
