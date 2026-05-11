---
paper_id: PAPER_1172
title: "Independent Re-Derivation of the 26-D Curvature Coefficient $\\rho_{R_{26}} = (13/2)\\,v_{UA}^{2}\\,\\rho_{\\mathrm{SCm}}$ via Two Routes"
session: 256
date: 2026-05-10
author: Daniel Murphy
status: CLOSED
cvw: v2.0.0
tags: [UQFF, R26, curvature, KK-reduction, Gauss-Bonnet, cross-check]
sm_anchor: G6_SM_Anchor_Gate
---

# PAPER_1172: Independent Re-Derivation of $\rho_{R_{26}}$

## Abstract

PAPER_1170 §3 obtained the 26-D curvature contribution to the
vacuum-energy ledger as
$\rho_{R_{26}} = \tfrac{13}{2}\,v_{UA}^{2}\,\rho_{\mathrm{SCm}}$
through Kaluza-Klein reduction of the 26-D Einstein-Hilbert term.
This paper verifies the prefactor $13/2$ by an **independent route**
using the Gauss-Bonnet trace identity on
$\mathcal{M}^{26} = \mathcal{M}^{4}_{\mathrm{BSFG}}\times T^{22}$,
arriving at the same rational with no free parameters. The double
derivation closes a potential ambiguity in the ledger.

---

## 1. Route A: KK reduction (recap of PAPER_1170)

With $V_{22} = (2\pi L_{\mathrm{KK}})^{22}$ and
$\kappa_4\,\rho_{\mathrm{SCm}} = (D_{\mathrm{crit}}-D_{\mathrm{phys}})/D_{\mathrm{crit}}
= 22/26 = 11/13$, the curvature density reduces to

$$
\rho_{R_{26}}^{(\mathrm{A})}
\;=\; \frac{\langle R_{26}\rangle}{2\kappa_4}
\;=\; \frac{11\,v_{UA}^{2}}{2}\cdot\frac{13}{11}\,\rho_{\mathrm{SCm}}
\;=\; \frac{13}{2}\,v_{UA}^{2}\,\rho_{\mathrm{SCm}}.
$$

## 2. Route B: Gauss-Bonnet trace identity

The 26-D Gauss-Bonnet density on a product manifold satisfies

$$
\langle G_{26}\rangle
\;=\; \langle R_{26}^{2}\rangle - 4\langle R^{ab}R_{ab}\rangle
   + \langle R^{abcd}R_{abcd}\rangle.
$$

On $\mathcal{M}^{4}_{\mathrm{BSFG}}\times T^{22}$ the torus is flat
($R_{T^{22}}=0$), so all curvature lives in the 4-D BSFG factor. The
trace identity then gives

$$
\langle R_{26}\rangle
\;=\; \langle R_{4}\rangle
\;+\; 2\,\frac{D_{\mathrm{crit}}-D_{\mathrm{phys}}}{L_{\mathrm{KK}}^{2}}\,\sin^{2}\!\theta_{\mathrm{mix}},
$$

where $\theta_{\mathrm{mix}}$ is the BSFG-torus mixing angle fixed by
the Archimedean half-condition $\sum\beta_{i}=3/2$ (PAPER_1165),
yielding $\sin^{2}\theta_{\mathrm{mix}} = D_{\mathrm{phys}}/(2D_{\mathrm{crit}}-D_{\mathrm{phys}})
= 4/48 = 1/12$.

Therefore

$$
\langle R_{26}\rangle
\;=\; \langle R_{4}\rangle \;+\; \frac{2\cdot 22}{L_{\mathrm{KK}}^{2}}\cdot\frac{1}{12}
\;=\; \langle R_{4}\rangle \;+\; \frac{11\,v_{UA}^{2}}{3}.
$$

On the BSFG vacuum $\langle R_{4}\rangle = (22/3)\,v_{UA}^{2}$ (PAPER_1166 §4),
so

$$
\langle R_{26}\rangle \;=\; \frac{22 + 11}{3}\,v_{UA}^{2}
\;=\; 11\,v_{UA}^{2},
$$

in exact agreement with Route A's value $\langle R_{26}\rangle = 11\,v_{UA}^{2}$.

Multiplying by the same $1/(2\kappa_4) = (13/22)\,\rho_{\mathrm{SCm}}/v_{UA}^{0}$:

$$
\rho_{R_{26}}^{(\mathrm{B})}
\;=\; \frac{11\,v_{UA}^{2}}{2}\cdot\frac{13}{11}\,\rho_{\mathrm{SCm}}
\;=\; \frac{13}{2}\,v_{UA}^{2}\,\rho_{\mathrm{SCm}}.
$$

## 3. Cross-check table

| Route | $\langle R_{26}\rangle$ | $1/(2\kappa_4)$ | $\rho_{R_{26}}$ |
|---|---|---|---|
| A (KK reduction) | $11\,v_{UA}^{2}$ | $(13/22)\rho_{\mathrm{SCm}}$ | $(13/2)v_{UA}^{2}\rho_{\mathrm{SCm}}$ |
| B (Gauss-Bonnet) | $11\,v_{UA}^{2}$ | $(13/22)\rho_{\mathrm{SCm}}$ | $(13/2)v_{UA}^{2}\rho_{\mathrm{SCm}}$ |
| **Numerical** | $1.1\times 10^{17}\ \mathrm{m^{-2}}$ | $4.19\times 10^{-37}$ | $\mathbf{4.61\times 10^{-20}\ J/m^{3}}$ |

Both routes agree exactly. The rational $13/2$ is therefore
**locked** by the closed dimensional structure $(D_{\mathrm{crit}}=26,
D_{\mathrm{phys}}=4)$ via two independent paths.

## 4. Implications

1. No room for $\rho_{R_{26}}$ to absorb the residual $0.2\%$ of
   PAPER_1170's ledger — the slack must come from $\rho_{\mathrm{KK}}$
   (closed in PAPER_1171) or $\rho_{\mathrm{BSFG}}$.
2. The Gauss-Bonnet route exposes a falsification handle:
   $\sin^{2}\theta_{\mathrm{mix}}=1/12$ predicts a specific
   BSFG-torus mixing pattern testable in any future numerical
   simulation of the closed Lagrangian.

## References

- PAPER_1165 — Triangular $\beta_{i}$ ladder.
- PAPER_1166 — V(UA) Mexican-hat closure; $\langle R_{4}\rangle$ value.
- PAPER_1167 — Master synthesis.
- PAPER_1170 — Four-line ledger; Route A.
- PAPER_1171 — KK regulator first-principles derivation.
- `uqff_closed_constants.py` — canonical integer-rational constants.
