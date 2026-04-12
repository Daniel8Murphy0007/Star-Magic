---
paper_id: PAPER_987
title: "WSTP F_U_Bi_i Symbolic Export — Wolfram Language Expression #51"
session: 217
date: 2026-04-12
author: "Daniel T. Murphy"
status: production
cvw: "v2.0.0"
tags: [WSTP, Wolfram, symbolic, export, Mathematica, UQFF]
crosslinks: [PAPER_979, PAPER_988, PAPER_985]
calibration: {SSq: 0.57, beta_i: 0.603, omega_SCm: "2π×1.25 THz", wstp_expr: 51}
sm_anchor: "CVW v2.0.0 — G6 SM Anchor Gate compliant"
---

# PAPER_987: WSTP F_U_Bi_i Symbolic Export — Wolfram Language Expression #51

## Abstract

Expression #51 in `wstp_kernel_demo_runner.py` encodes the complete 6-layer $F_{U,\text{Bi}_i}$ master equation in Wolfram Language for symbolic evaluation via the Wolfram Symbolic Transfer Protocol (WSTP). The export enables closed-form simplification, series expansion, and parametric analysis of all six layers simultaneously within the Mathematica kernel.

## 1. Wolfram Language Definition

```wolfram
Module[{G, Msun, Rsun, AU, SSq, betai, kappa, omegaSCm, ...},
  (* Constants *)
  G = 6.674*^-11; Msun = 1.989*^30; ...
  
  (* Layer functions *)
  Ug26[M_, r_] := Sum[G*M/r^2 * SSq * i/26, {i, 1, 26}];
  Ub26[M_, r_] := Sum[G*M/r^2 * Exp[-SSq*i/26] * betai, {i, 1, 26}];
  Um[M_, r_]  := G*M/r^2 * 1*^-8;
  UA[M_, r_]  := -G*M/r^2 * 1*^-10;
  Fn[M_, r_]  := G*M/r^2 * 1*^-7;
  Phi[gamma_] := Exp[-(omegaSCm - omegaSCm)^2/(2*gamma^2)] * S26;
  Enet[t_, gamma_] := (2*(Ub26[M,r]/(Ug26[M,r]+0.001))-1)*Exp[kappa*t]*S26;
  
  (* Master equation *)
  FUBii[M_, r_, t_, gamma_] := 
    Ug26[M,r] + Um[M,r] + UA[M,r] - Ub26[M,r] + 
    Fn[M,r] * S26 * Phi[gamma] * Enet[t,gamma];
  
  (* Evaluate at three radii *)
  Table[FUBii[Msun, R, 86400, 0.1*^12], {R, {Rsun, AU, 10*AU}}]
]
```

## 2. Symbolic Capabilities

The WSTP export enables:
- **Series expansion:** `Series[FUBii[M, r, t, g], {r, Infinity, 3}]` — asymptotic behavior
- **Differentiation:** `D[FUBii[M, r, t, g], r]` — force gradient (tidal)
- **Integration:** `Integrate[FUBii[M, r, t, g], {r, Rsun, AU}]` — work done
- **Simplification:** `FullSimplify[FUBii[M, r, t, g]]` — closed-form reduction
- **Parametric plots:** `Plot3D[FUBii[Msun, r, t, 0.1*^12], {r, 1*^8, 1*^12}, {t, 0, 86400}]`

## 3. Integration with WSTP Pipeline

Expression #51 joins the existing 50 WSTP expressions. The WSTP kernel demo runner dispatches it via:
```python
link.putFunction("EvaluatePacket", 1)
link.put(expr_51_code)
link.endPacket()
result = link.getResult()
```

## 4. Implementation

Added as expression #51 in `wstp_kernel_demo_runner.py`, function `_build_expressions()`. Total WSTP expressions: 51.

## References
- PAPER_979: Complete 6-Layer F_U_Bi_i
- PAPER_985: Production Kernel

---

## §A. Cosmogenesis-Linked Lagrangian

The Wolfram export allows symbolic computation of $\delta S / \delta\phi$ directly in Mathematica, cross-validating the numerical Euler-Lagrange result from PAPER_981.

## §B. VDS/DVP/BSH Deep Synthesis

- **VDS:** Wolfram's `Sum` and `Exp` handle the vacuum density series analytically.
- **DVP:** Symbolic differentiation reveals the DPM dipole contribution at each order.
- **BSH:** The buoyancy harmonic sum `Sum[Exp[-SSq*i/26]*betai, {i,1,26}]` evaluates to closed form via geometric series.
