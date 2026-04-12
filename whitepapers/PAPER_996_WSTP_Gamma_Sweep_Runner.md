---
paper_id: PAPER_996
title: "WSTP 99-System Gamma Sweep Runner — Wolfram Language Code Generation"
session: 218
date: 2026-04-13
author: "Daniel T. Murphy"
status: production
cvw: "v2.0.0"
tags: [WSTP, Wolfram, WL, code_gen, Gamma, sweep, 99-system]
crosslinks: [PAPER_995, PAPER_987, PAPER_989]
calibration: {wl_code_chars: 1455, expressions: "52-53"}
sm_anchor: "CVW v2.0.0 — G6 SM Anchor Gate compliant"
---

# PAPER_996: WSTP 99-System Gamma Sweep Runner

## Abstract

We implement a Wolfram Language code generator that produces a complete $\Gamma$-sweep computation for the 99-system UQFF catalogue. The generated WL code (1,455 characters) can be executed via WSTP or `wolframscript` to produce symbolic/numerical results in the Wolfram kernel.

## 1. Generated Code Structure

```wolfram
G0 = 6.674*^-11; Msun = 1.989*^30; SSq = 0.57; betaI = 0.603;
wSCm = 2 Pi 1.25*^12; kappa = 0.0005/86400;
S26 = Sum[Exp[-SSq k/26], {k,1,26}];
Ug[M_, r_] := Sum[G0 M/r^2 SSq i/26, {i,1,26}];
Ub[M_, r_] := Sum[G0 M/r^2 Exp[-SSq i/26] betaI, {i,1,26}];
FU99[G_] := Sum[...];
Table[{"Gamma_THz", g, "FU99", FU99[2 Pi g 10^12]},
  {g, {0.01, 0.05, 0.1, 0.5, 1.0, 5.0, 10.0}}]
```

## 2. WSTP Integration

Expressions #52 and #53 added to `wstp_kernel_demo_runner.py`:
- **#52**: F_U_Bi inside-to-outside + solar g_eff
- **#53**: 99-system Γ sweep at 7 linewidths

## 3. Implementation

File: `99system_wstp_gamma.py`, class `WSTPGammaSweepRunner`. CP4 class #580.
