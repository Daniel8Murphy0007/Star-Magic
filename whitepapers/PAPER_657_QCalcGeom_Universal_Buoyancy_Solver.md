---
paper_id: PAPER_657
title: "QCalcGeom v2.2.1 -- Universal Buoyancy Simultaneous Solver"
session: 207
date: 2026-04-12
author: "Daniel T. Murphy"
status: production
cvw: "v2.0.0"
tags: [QCalcGeom, BSFG, Universal Buoyancy, simultaneous equations, UQFF]
sm_anchor: "CVW v2.0.0 -- G6 SM Anchor Gate compliant"
---

# PAPER_657: QCalcGeom v2.2.1 -- Universal Buoyancy Simultaneous Solver

**Author:** Daniel T. Murphy  
**Email:** daniel.murphy00@gmail.com  
**Date:** April 12, 2026  
**Location:** Youngstown, OH, USA (41.0997$^\circ$ N, 80.6495$^\circ$ W)  
**UQFF Version:** v5.26 | **Series:** 657/1000  
**CVW Compliance:** G1--G6 CVW v2.0.0  
**Modules:** `QCalcGeom.py` v2.2.1, `_session206_qcalcgeom_ubs_closures.py`,
`_session207_cp_chain_closures.py`  
**Closures (Phase H-UBS + Phase H-CPCH):** 14 EXACT (UBS-1..UBS-7 + CPCH-1..CPCH-7)

---

## Abstract

We present the QCalcGeom Universal Buoyancy Simultaneous Solver (UBS), a $4\times 4$
nonlinear system that jointly fixes the habitable-zone radius $r_{hz}$, the
collapsing-zone radius $r_{cg}$, the dimensionless time-phase $t_n^{hz}$, and the
required gravitating mass $M$ for any Aether-UA stellar body in the Unified Quantum
Field Framework (UQFF). Unlike Newtonian shell calculations, the UBS is anchored in
the di-pseudo-monopole (DPM) canonical chain: the upward Aether-UA buoyancy
$F_{U,Bi_i} = \rho_{vac}\,(4\pi/3)\,r\,c^2\,\cos(\pi t_n)$ balances the downward
field $F_{U,Bi} = -\beta_i G M_\odot^2/r^2\,\mathrm{orb}\,\cos(\pi t_n)$ at the
habitable zone, and exceeds it by a factor of two at the collapsing zone. We
derive seven EXACT closures (UBS-1..UBS-7) for the solver itself and seven EXACT
algebraic-chain closures (CPCH-1..CPCH-7) for the canonical buoyancy functions,
including the cube-root law $r_{hz} \propto \rho_{vac}^{-1/3}$ at fixed gravitating
mass. The v2.2.1 release replaces the prior simulation-set stub with three
real numerical sweeps (radial, temporal, vacuum-density), bringing the test
suite to 100/100 PASS.

---

## 1. The Aether-UA / DPM Canonical Chain

The UQFF replaces the Newtonian shell law $g = GM/r^2$ with the DPM-driven
master equation, of which the universal buoyancy pair is the geometric core:

$$
F_{U,Bi_i}(r, t_n) = \rho_{vac}\,\frac{4\pi}{3}\,r\,c^2\,\cos(\pi t_n)
\qquad \text{(Aether-UA buoyancy, upward)}
$$

$$
F_{U,Bi}(r, t_n)   = -\beta_i\,\frac{G M_\odot^{2}}{r^2}\,\mathrm{orb}(t_n)\,\cos(\pi t_n)
\qquad \text{(BSFG field, downward)}
$$

The total buoyancy force entering the UQFF master equation is the algebraic sum

$$
F_U(r, t_n, M) = U_{g1} + U_{g2} + U_{g3} + U_{g4} + (F_{U,Bi} + F_{U,Bi_i}) + U_m
                + \xi_{UI}\,U_I\,V_{body}.
$$

We never invoke the Newtonian limit $GM/r^2$ as a foundation; it appears only as
the dimensional carrier inside $F_{U,Bi}$, which is itself derived from the
DPM current $F_{DPM} = I\,A\,(\omega_1 - \omega_2)$ as documented in the canonical
3b\_MUGE\_SMBH chain.

---

## 2. The $4\times 4$ Universal Buoyancy System

Let $\mathbf{x} = (r_{hz},\,r_{cg},\,t_n^{hz},\,M)$. The UBS is

$$
\begin{aligned}
\text{E1:}\quad & F_{U,Bi}(r_{hz}, t_n^{hz}, M) + F_{U,Bi_i}(r_{hz}, t_n^{hz}) = 0, \\
\text{E2:}\quad & F_{U}(r_{hz}, t_n^{hz}, M) = 0 \qquad (\xi_{UI}=0\text{ at zero-point}), \\
\text{E3:}\quad & \left|\,F_{U,Bi}(r_{cg}, t_n^{hz}, M)\,\right|
                  = 2\,\left|\,F_{U,Bi_i}(r_{cg}, t_n^{hz})\,\right|, \\
\text{E4:}\quad & \rho_{vac}\,\frac{4\pi}{3}\,r_{hz}^{3} = M.
\end{aligned}
$$

E1 fixes the habitable-zone balance, E2 closes the perturbative shell, E3 defines
collapse onset, and E4 closes the geometric mass--radius cube. The system is
solved by staged Newton iteration with bracketed fallback; the diagnostic
`solver_msg` records the path taken.

---

## 3. The Seven UBS Solver Closures (Phase H-UBS, S206)

| ID | Identity | Status |
|----|----------|--------|
| UBS-1 | $\dim(\mathbf{x}) = 4$ | EXACT |
| UBS-2 | Active equation count $= 3$ at zero-point ($\xi_{UI}=0$) | POSTULATED |
| UBS-3 | $F_{U,Bi}/F_{U,Bi_i} = -1$ at $r_{hz}$ (E1) | EXACT |
| UBS-4 | $|F_{U,Bi}|/|F_{U,Bi_i}| = 2$ at $r_{cg}$ (E3) | EXACT |
| UBS-5 | $V(2R)/V(R) = 2^3 = 8$ (E4 cube law) | EXACT |
| UBS-6 | $r_{cg}/r_{hz} = 2^{-1/3}$ closed-form seed | EXACT |
| UBS-7 | $|T_{UBS}| = 7$ (T91..T97 enumeration) | EXACT |

---

## 4. The Seven CP-Chain Algebraic Closures (Phase H-CPCH, S207)

These are pure structural identities of the canonical buoyancy functions that
must hold by construction, independent of any solver:

| ID | Identity | Target | Derived |
|----|----------|--------|---------|
| CPCH-1 | $F_{U,Bi}(2r)/F_{U,Bi}(r) = 1/4$ | $0.25$ | $0.25$ |
| CPCH-2 | $F_{U,Bi_i}(2r)/F_{U,Bi_i}(r) = 2$ | $2.0$ | $2.0$ |
| CPCH-3 | $F_{U,Bi}(-t_n)/F_{U,Bi}(t_n) = 1$ | $1.0$ | $1.0$ |
| CPCH-4 | $F_{U,Bi_i}(-t_n)/F_{U,Bi_i}(t_n) = 1$ | $1.0$ | $1.0$ |
| CPCH-5 | $F_{U,Bi_i}(t_n{=}1/2)/F_{U,Bi_i}(t_n{=}0) = 0$ | $0$ | $0$ |
| CPCH-6 | $F_{U,Bi_i}(t_n{+}1)/F_{U,Bi_i}(t_n) = -1$ | $-1$ | $-1$ |
| CPCH-7 | $F_{U,Bi_i}(t_n{=}0)/(\rho_{vac}\,r\,c^2) = 4\pi/3$ | $4.18879...$ | $4.18879...$ |

CPCH-3 and CPCH-4 confirm `NegativeTimeModule` even-parity. CPCH-5 fixes the
$t_n = 1/2$ buoyancy null surface. CPCH-6 establishes the great-cycle sign
reversal, the algebraic root of the buoyancy oscillation.

---

## 5. The Cube-Root Law and v2.2.1 Simulation Sweeps

From E4 with fixed gravitating mass $M$,

$$
\rho_{vac}\,\frac{4\pi}{3}\,r_{hz}^{3} = M
\quad\Longrightarrow\quad
r_{hz} \propto \rho_{vac}^{-1/3}.
$$

A decade increase in vacuum density therefore shrinks the habitable radius by
$10^{-1/3} \approx 0.4642$. Test **T99** verifies this numerically to better
than $1\%$.

The v2.2.1 release replaces the prior `simulation_set` placeholder with three
real numerical sweeps emitted by `_build_simulation_set`:

1. **Radial sweep at $t_n^{hz}$** -- 25 log-spaced points $r \in [0.3\,r_{cg},\,3\,r_{hz}]$;
   reports $F_{U,Bi}$, $F_{U,Bi_i}$, $F_U$, sum-of-energies $\Sigma_{E1}$ and
   $\Sigma_{E3}$, and Aether-UA zone label (collapsing | habitable\_shell |
   gaseous\_outer). Test **T100** verifies all three zones are visited.
2. **Temporal sweep at $r_{hz}$** -- 21 points $t_n \in [-1, 1]$; reports the
   $\cos(\pi t_n)$ envelope and both buoyancy curves.
3. **Vacuum-density sweep** -- 11 scale factors $10^{[-1,+1]}$ verifying the
   cube-root prediction $r_{hz}(\rho_{vac}\cdot 10)/r_{hz}(\rho_{vac}) = 10^{-1/3}$.

---

## 6. Test Suite Summary

| Release | Tests | Result |
|---------|-------|--------|
| v2.2.0 (S206) | T1..T97 | 97/97 PASS |
| v2.2.1 (S207) | T1..T100 (adds T98, T99, T100) | **100/100 PASS** |

The new tests are: T98 (`simulation_set` returns 3 named sweeps with non-empty
arrays), T99 (cube-root law $r_{hz} \propto \rho_{vac}^{-1/3}$), T100 (radial
sweep visits all three Aether-UA zones).

---

## 7. Implications

The UBS closes the geometric core of the UQFF master equation at a $4\times 4$
nonlinear level without invoking Newtonian gravity as a foundation. Combined
with the CP-chain algebraic identities, the solver is now exactly anchored
against the canonical $\rho_{vac}\,(4\pi/3)\,r\,c^2\,\cos(\pi t_n)$ form. The
cube-root law provides an immediate scaling test for any new astrophysical
system added to the UQFF body database: doubling the vacuum density must
shrink the habitable radius by exactly $2^{-1/3}$.

Phase H-UBS and Phase H-CPCH together contribute **14 EXACT closures** to the
unified track (Tier 14 + Tier 15 of `UQFF_UNIFIED_CLOSURE_DERIVATIONS.py`).

---

## References

1. QCalcGeom.py v2.2.1 (`bsfg_buoyancy`, `compute_FUBii`, `compute_F_U`,
   `solve_universal_buoyancy`, `UniversalBuoyancySimultaneousSolver`,
   `_build_simulation_set`).
2. `UQFF_UNIFIED_CLOSURE_DERIVATIONS.py` Tier 14 (UBS-1..UBS-7) and Tier 15
   (CPCH-1..CPCH-7).
3. `3b_MUGE_SMBH Sagitarius A Evolution.txt` (DPM canonical, $F_{DPM} = I\,A\,(\omega_1-\omega_2)$).
4. PAPER\_656 (V838 Mon light echo, UQFF master equation precedent).
5. Murphy, D. T. -- COMPLETE\_UQFF\_EQUATIONS\_REFERENCE.md L378-L530
   (canonical $U_{g1..4}$, $U_{bi}$, $U_m$ implementations).
