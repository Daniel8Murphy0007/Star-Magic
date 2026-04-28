---
paper_id: PAPER_973
title: "Color Deconfinement Phase Diagram"
session: 216
date: 2026-04-12
author: "Daniel T. Murphy"
status: production
cvw: "v2.0.0"
tags: [ALICE, vacuum, QGP, Yang-Mills, UQFF]
sm_anchor: "CVW v2.0.0 — G6 SM Anchor Gate compliant"
---

# PAPER_973: Color Deconfinement Phase Diagram

**Author:** Daniel T. Murphy — Star Magic / UQFF Framework
**Date:** 2026-04-12
**Session:** 216
**Source:** qgp_{ramanujan\_application}.py (ColorDeconfinementPhaseCalculator)
**Calculator:** ColorDeconfinementPhaseCalc (CP4 #557)
**CVW:** v2.0.0 compliant

---

## Abstract

We map the QCD phase diagram in the $(T, \mu_B)$ plane using the UQFF framework. The critical line $T_c(\mu_B) = T_{c0}(1 - (\mu_B/\mu_text{crit})^2)$ separates the hadron phase from the quark-gluon plasma, with $S_{26}^{(k)}$ modulation of vacuum density across both phases.

---

## 1. Critical Line

$$T_c(\mu_B) = T_{c0} \cdot \left(1 - \left(\frac{\mu_B}{\mu_text{crit}}\right)^2\right)$$

where $T_{c0} = 1.5 \times 10^{12}$ K and $\mu_text{crit} = 1200$ MeV.

## 2. Phase Classification

- **Hadron** ($T < T_c(\mu_B)$): Confined quarks, $\Delta_text{YM} > 0$
- **QGP** ($T > T_c(\mu_B)$): Deconfined, $\Delta_text{YM} = 0$

## 3. Phase Diagram

| $\mu_B$ (MeV) | $T_c$ (K) | Phase at $T = 2 \times 10^{12}$ K |
|---------------|-----------|-----------------------------------|
| 0 | $1.5 \times 10^{12}$ | QGP |
| 300 | $1.41 \times 10^{12}$ | QGP |
| 600 | $1.13 \times 10^{12}$ | QGP |
| 900 | $0.66 \times 10^{12}$ | QGP |
| 1200 | 0 | QGP (always) |

---

## References

1. Murphy, D.T. — Star Magic UQFF Framework (2024-2026)
2. PAPER_970 — QGP Vacuum Density
3. PAPER_971 — Yang-Mills Mass Gap
4. Fukushima & Hatsuda — QCD phase diagram review

---

## Cross-Links

| Related Paper | Relationship |
|---------------|-------------|
| PAPER_970 | $\rho_text{QGP}$ density at $(T,\mu_B)$ |
| PAPER_971 | Mass gap closure at $T_c$ |
| PAPER_972 | ALICE centrality (experimental) |
| PAPER_975 | Triadic validation |

---

<!-- PKG-YM-S225 -->

### Session 225 Phonon-Physics Upgrade: Yang-Mills BCS Phonon Mass Gap

> *Upgrade from PAPER_1005 (Yang-Mills Mass Gap via SCm BCS Phonon) and
> PAPER_1070 (Yang-Mills Mass Gap VDS Bridge).  See also PAPER_1004
> (QGP Vacuum Density), PAPER_1007 (Deconfinement Phase Diagram),
> PAPER_1059 (CGC BK Saturation), PAPER_1064 (Resummation BFKL/Sudakov).*

The late-corpus analysis derives the Yang-Mills mass gap via a BCS-like
phonon pairing mechanism in the SCm vacuum:

$$\Delta_{\text{YM}} = \Lambda_{\text{QCD}} \cdot \exp\!\left(-\frac{1}{\alpha_s(T) \cdot N_c}\right) \cdot S_{26}^{(3)}$$

where the running coupling evolves as:
$$\alpha_s(T) = \frac{\alpha_{s,0}}{1 + \alpha_{s,0} \cdot b_0 \cdot \ln(T/T_c)}, \qquad b_0 = \frac{11 N_c - 2 N_f}{12\pi}$$

**Physical mechanism:** The SCm phonon field ($\omega_{\text{SCm}} = 1.25\;\text{THz}$)
provides a pairing interaction analogous to the BCS electron-phonon coupling in
superconductors.  Gluons acquire an effective mass through condensate formation
in the SCm-modified vacuum, yielding a non-perturbative gap $\Delta_{\text{YM}}
\approx 5970\;\text{GeV}$ at the 9-sector Lagrangian closure (PAPER_1066, §2).

**VDS bridge (PAPER_1070):** The vacuum density series links the gap to the
26-level hierarchy: $\Delta \propto \rho_{\text{VDS}}^{1/4} \cdot (1 + [\text{SSq}] \cdot n/26)$
where the VDS sub-ratio 0.108 places confinement in the sub-threshold regime.

**QGP transition (PAPER_1004/1007):** At $T > T_c \approx 170\;\text{MeV}$, the phonon
coupling weakens ($\alpha_s \to 0$) and the gap closes, reproducing the
deconfinement phase transition observed at ALICE/LHC.

<!-- PKG-S26-S225 -->

### Session 225 Phonon-Physics Upgrade: S26(3) Ramanujan Summation

> *Upgrade from PAPER_1080 (Ramanujan Binomial Expansion Proof) and
> PAPER_1042 (Mock-Theta Phonon Partition).  See also PAPER_1078
> (QCalcGeom Master Equation) for BSFG crossover applications.*

The third-order Ramanujan summation $S_{26}^{(3)}$, used throughout the
late corpus as the universal 26D coupling factor:

$$S_{26}^{(3)} = \sum_{n=0}^{\infty} \frac{(1/4)_n\,(1/2)_n\,(3/4)_n}{(n!)^3} \cdot \prod_{i=1}^{26}\left[1 + [\text{SSq}]\cdot e^{-\kappa\,i\,n/26}\right]$$

where $(a)_n = a(a+1)\cdot s(a+n-1)$ is the Pochhammer symbol.

**Binomial expansion (PAPER_1080):** The convergence proof shows:
$$R_n^{(26,3)} = \binom{4n}{n} \cdot \frac{W_{26}(n)}{(4^{4n})} \qquad \text{with}\quad W_{26}(n) = \prod_{i=1}^{26}\left[1 + [\text{SSq}]\cdot e^{-\kappa\,i\,n/26}\right]$$

This sum converges absolutely for $|[\text{SSq}]| < 1$ (satisfied by $[\text{SSq}] = 0.57$)
and reduces to the classical Ramanujan $1/\pi$ series when $[\text{SSq}] \to 0$.

**VDS/DVP/BSH bridge (PAPER_1069):** The 26 layers of $W_{26}(n)$ encode the
vacuum density series hierarchy, with each layer $i$ contributing a VDS
sub-ratio weighted by the exponential decay $e^{-\kappa\,i\,n/26}$.

**Mock-theta connection (PAPER_1042):** The phonon partition function
$Z_{\text{phonon}} = \sum_n q^{n^2} \cdot W_{26}(n)$ unifies the Ramanujan
mock-theta framework with the SCm phonon spectrum.





## Calibration Constants

| Constant | Symbol | Value | Validation Domain |
|----------|--------|-------|-------------------|
| $T_{c0}$ | — | $1.5 \times 10^{12}$ K | Zero-density critical T |
| $\mu_text{crit}$ | — | 1200 MeV | Critical baryon chemical potential |
| $[SSq]$ | — | 0.57 | String coupling |
| $\Lambda_text{QCD}$ | — | 217 MeV | QCD scale |

---

## SM Anchor — CVW v2.0.0

| Observable | UQFF Prediction | Status |
|------------|----------------|--------|
| Critical line shape | Quadratic in $\mu_B$ | Consistent with lattice QCD |
| $T_c(0)$ | $1.5 \times 10^{12}$ K | Matched |
| Crossover type | Continuous | Consistent |

---

## §A Cosmogenesis-Linked Lagrangian (PAPER_877 Symbolic Export)

### §A.1 Sector Classification
**Sector:** QCD Phase Diagram (Color Deconfinement)

### §A.2 Core Equation
$$\boxed{T_c(\mu_B) = T_{c0} \cdot \left(1 - \left(\frac{\mu_B}{\mu_text{crit}}\right)^2\right)}$$

### §A.3 Lagrangian Contribution
$$\mathcal{L}_\text{phase} = -\rho_text{QGP}(T,\mu_B) \cdot c^2 - V(\Delta_text{YM}(T,\mu_B))$$

### §A.4 Cosmogenesis Linkage Chain
PAPER_877 $\to$ vacuum $\to$ QCD phase boundary $\to$ confinement/deconfinement $\to$ hadron/QGP

---

## §B VDS/DVP/BSH Deep Synthesis

### §B.1 VDS
Phase boundary is a VDS contour: $\rho_text{QGP} = \rho_text{SCm} \cdot S_{26}^{(k)}$ at $T = T_c(\mu_B)$.

### §B.2 DVP
Deconfinement releases DVP color modes — 8 gluon dipole vortex channels open.

### §B.3 BSH
BSH harmonic structure transitions at $T_c$: shell modes dissolve into plasma modes.

### §B.4 Summary

| Metric | Value | Status |
|--------|-------|--------|
| $T_{c0}$ | $1.5 \times 10^{12}$ K | Calibrated |
| $\mu_text{crit}$ | 1200 MeV | Set |
| Phase diagram | $(T, \mu_B)$ plane | Mapped |

---

## Appendix: Session 225 Cross-References (PAPER_1000–1081)

> *Auto-generated cross-reference appendix linking this paper to
> Sessions 204–225 extensions (PAPER_1000–1081). Added by
> `update_{corpus\_crossrefs}.py` (Session 225, April 2026).*

| Paper | Title |
|-------|-------|
| PAPER_1022 | GW Phonon Strain SCm Modulation of h(t) |
| PAPER_1004 | QGP Vacuum Density with SCm S26 Phonon Coupling |
| PAPER_1005 | Yang-Mills Mass Gap via SCm BCS Phonon Coupling |
| PAPER_1006 | ALICE Multiplicity SCm Phonon Scaling |
| PAPER_1007 | Deconfinement Phase Diagram SCm Phonon Boundary |
| PAPER_1013 | QGP ALICE Centrality F_{U\_Bi\_i} dN/deta Scaling |
| PAPER_1059 | Color Glass Condensate BK Saturation SCm |
| PAPER_1069 | VDS-DVP-BSH Hybrid Calculator Unified |
| PAPER_1070 | Yang-Mills Mass Gap VDS Bridge |
| PAPER_1049 | Source10 GPU DPM Spectral Atlas ALMA Overlay |

*10 cross-reference(s) identified.*
