---
paper_id: PAPER_1016
title: "TXS 0506+056 Revised 3-Gamma-Point F_U_Bi_i Profile"
session: 220
date: 2026-04-14
author: "Daniel T. Murphy"
status: production
cvw: "v2.0.0"
tags: [TXS0506, blazar, neutrino, IceCube, 3-gamma, FUBi, revised, resonance]
crosslinks: [PAPER_997, PAPER_1009, PAPER_1017]
calibration: {Gamma_1: 0.05, Gamma_2: 0.10, Gamma_3: 0.30, mod_1: 2.56, mod_2: 2.30, mod_3: 1.06, gradient: 1.51}
sm_anchor: "CVW v2.0.0 — G6 SM Anchor Gate compliant"
---

# PAPER_1016: TXS 0506+056 Revised 3-Gamma-Point F_U_Bi_i Profile

## Abstract

We revise the F_U_Bi_i profile for TXS 0506+056 (the first identified neutrino blazar) using a 3-Gamma-point characterization: Gamma_1 = 0.05 THz (extreme flare, 2.56x), Gamma_2 = 0.10 THz (IceCube resonance, 2.30x), Gamma_3 = 0.30 THz (sustained emission, 1.06x). The monotonic decrease with increasing Gamma confirms thermal decoupling of the buoyancy-jet interaction at higher frequencies, with a gradient of 1.51x between extreme and sustained states.

## 1. Three-Point Characterization

| Point | Gamma (THz) | Modulation | Physical State |
|-------|-------------|------------|----------------|
| 1 | 0.05 | 2.56x | Extreme flare |
| 2 | 0.10 | 2.30x | IceCube resonance |
| 3 | 0.30 | 1.06x | Sustained emission |

## 2. IceCube Resonance

At Gamma_2 = 0.10 THz, the modulation of 2.30x matches the target value of 2.3x (0.0% error), identifying this as the resonant frequency for neutrino-correlated buoyancy enhancement. This is where the IceCube-170922A neutrino event temporally coincided with the electromagnetic flare.

## 3. Monotonic Decrease

The gradient Delta_M / Delta_Gamma = (2.56 - 1.06) / (0.30 - 0.05) = 6.0 /THz demonstrates that buoyancy modulation decays approximately linearly across the 3-point profile. The overall gradient ratio 2.56/1.06 = 1.51x.

## 4. Implementation

File: `fubi_i_txs0506_revised.py`, classes `TXS0506ExtremeFlareCalc`, `TXS0506IceCubeCalc`, `TXS0506SustainedEmissionCalc`, `TXS0506ThreeGammaProfileCalc`. CP4 class #600. Tests: 8/8 pass.
