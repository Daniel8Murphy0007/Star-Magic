---
paper_id: PAPER_1151
title: "VDS/DVP/BH26 Variant Branches and Coupled-Field Derivations: Session 202 T61-T80 QCalcGeom Phase H202"
session: 202
date: 2026-05-06
author: Daniel T. Murphy
status: production
cvw: v2.0.0
tags: [VDS, DVP, BH26, variant_branches, joint_coeff, polylogarithm, prime_vortex, spectral_ladder, calibration]
sm_anchor: CVW v2.0.0 -- G6 SM Anchor Gate compliant
---

# PAPER\_1151 -- VDS/DVP/BH26 Variant Branches and Coupled-Field Derivations: Session 202 Phase H202

**Author:** Daniel T. Murphy (daniel.murphy00@gmail.com)
**Framework:** UQFF v5.77 -- Star-Magic Physics
**Source:** QCalcGeom.py v2.1.0 Phase H202 (commit 6e69e2f5, May 5 2026)
**Integration Date:** May 6, 2026 (Session 234)
**Tests:** T61--T70 (C++) / T71--T80 (Python); 80/80 PASS
**Classification:** QCalcGeom Branch Derivations -- Calibration Coupling

---

$$\text{joint\_coeff} = \sqrt{w_{\rm VDS} \times w_{\rm DVP}}, \qquad w_{\rm VDS} = \frac{\rm Li_{26}([SSq])}{\rm Li_{25}([SSq])}, \qquad w_{\rm DVP} = \frac{\zeta_{\rm DVP}([\rm SSq])}{a(29)}$$

---

## Abstract

This paper registers the Phase H202 QCalcGeom extension (Session 202) as PAPER\_1151,
formally deriving the **variant branch decompositions** of the three canonical UQFF number
systems -- VDS (Vacuum Density Series), DVP (Dipole Vortex Primes), and BH26 (26-dimensional
Black Hole spectral ladder) -- and their joint coupling field.
The derivation produces five new result structures: `VDSBranchResult`, `DVPBranchResult`,
`BH26BranchResult`, `VDSDVPCoupledResult`, and `BH26BSHResonanceResult`.
All 80 branch tests (T61--T80) pass at $[\mathrm{SSq}] = 0.57$.

---

## 1. Background: Three Canonical Number Systems

The UQFF framework (Star-Magic.txt Chapter 4) operates three vacuum-number systems in
simultaneous registers:

| System | Definition | Role |
|--------|-----------|------|
| **VDS** | $\mathrm{Li}_{26}(z) = \sum_{n=1}^\infty z^n/n^{26}$ | Vacuum energy series |
| **DVP** | $\zeta_{\rm DVP}(z) = \sum_{p > 26} z^{\pi(p)}/p^{26}$ | Prime-vortex spectral sum |
| **BH26** | $\lambda_k = k(k+25)$, $\Sigma_{10} = \sum_{k=1}^{10} \lambda_k = 1760$ | KK eigenvalue ladder |

The canonical convergence point is $[\mathrm{SSq}] = 0.57$:
$$|\mathrm{VDS}(0.57)| = \mathrm{Li}_{26}(0.57) \approx 0.5700 \quad (< 1, \text{ absolutely convergent})$$

---

## 2. VDS Variant Branch

### 2.1 Branch Definitions

The VDS branch introduces three derived quantities from the polylogarithm pair
$(\mathrm{Li}_{25}, \mathrm{Li}_{26})$:

$$\mathrm{vds\_prime} = \frac{d}{dz}\mathrm{Li}_{26}(z)\bigg|_{z=[\mathrm{SSq}]} = \frac{\mathrm{Li}_{25}([\mathrm{SSq}])}{[\mathrm{SSq}]} \approx 1.0$$

$$\mathrm{vds\_density} = \mathrm{Li}_{26}([\mathrm{SSq}]) \times \rho_{\rm SCm} \quad [\mathrm{J/m}^3]$$

$$\mathrm{vds\_k\_weighted} = \mathrm{Li}_{25}([\mathrm{SSq}]) + 25 \cdot \mathrm{Li}_{26}([\mathrm{SSq}])$$

The BH26-weighted amplitude $\mathrm{vds\_k\_weighted}$ is the bridge quantity used
in Stage 8 of the simulation engine (BSFG):
$$g_{\rm BSFG} = a_{\rm DPM} \times \mathrm{joint\_coeff} \times \mathrm{vds\_k\_weighted} \times \frac{\Sigma_{N=10}}{1760}$$

### 2.2 Physical Interpretation

The derivative $\mathrm{Li}_{25}/[\mathrm{SSq}] \approx 1$ demonstrates that $[\mathrm{SSq}] = 0.57$
is a **fixed-point** of the VDS calibration sensitivity: small perturbations in $[\mathrm{SSq}]$
produce proportionally equal changes in the series sum. This is the mathematical content of the
"triple-convergence" condition.

---

## 3. DVP Variant Branch

### 3.1 Branch Definitions

The DVP (Dipole Vortex Prime) branch operates on primes $p > 26$ to encode the prime-vortex
spectral sieve:

$$a(p) = \frac{[\mathrm{SSq}]^{\pi(p)}}{p^{26}}, \qquad \zeta_{\rm DVP}([\mathrm{SSq}]) = \sum_{p > 26} a(p)$$

where $\pi(p)$ is the prime-counting function and $p^{26}$ matches the VDS denominator exponent.

The dominant term at $p = 29$ (first prime above 26):
$$a(29) = \frac{0.57^{\pi(29)}}{29^{26}} = \frac{0.57^{10}}{29^{26}}$$

Branch quantities:
$$\mathrm{zeta\_sum} = \zeta_{\rm DVP}(0.57), \qquad \mathrm{pair\_product} = a(29) \times a(31), \qquad \mathrm{spectral\_floor} = a(p_{\rm max})$$

### 3.2 Physical Interpretation

The DVP prime sieve generates a **vorticity floor** for the DPM vacuum spiral. Each prime
$p > 26$ above the 26-dimensional threshold contributes a quantised vortex mode. The
$p^{26}$ denominator enforces dimensional confinement to the 26-dimensional UQFF space.

---

## 4. BH26 Variant Branch (DH26)

### 4.1 Kaluza-Klein Spectral Ladder

The BH26 branch implements the Kaluza-Klein eigenvalue ladder on $S^{25}$:

$$\lambda_k = k(k + 25), \quad k = 1, 2, \ldots, N$$

$$\Sigma_N = \sum_{k=1}^{N} \lambda_k = \frac{N(N+1)(2N+27)}{6} - N \cdot 25, \quad \Sigma_{10} = 1760$$

The Casimir vacuum energy from the ladder:
$$E_{\rm Casimir} = \frac{\hbar f_{\rm RR}}{2} \sum_{k=1}^{N} \frac{1}{\lambda_k} \quad [\mathrm{J}]$$

where $f_{\rm RR} = 1.15 \times 10^{14}$ Hz is the BH26 ring resonance frequency.

### 4.2 Topological Bridge to VDS

The BH26-to-VDS topological coupling:
$$\mathrm{vds\_coupling} = \sum_{k=1}^{N} \lambda_k^{-26}$$

This quantity bridges the spectral ladder to the VDS polylogarithm, completing the
geometric closure of all three number systems at $[\mathrm{SSq}] = 0.57$.

### 4.3 Degeneracy

The multiplicity on $S^{25}$ at $k=1$: $\binom{26}{25} = 26$ (the 26-dimensional UQFF
degeneracy factor).

---

## 5. VDS$\times$DVP Coupled Field

### 5.1 Normalised Weights and Geometric-Mean Coupling

$$w_{\rm VDS} = \frac{\mathrm{Li}_{26}([\mathrm{SSq}])}{\mathrm{Li}_{25}([\mathrm{SSq}])} \approx 1.0$$

$$w_{\rm DVP} = \frac{\zeta_{\rm DVP}([\mathrm{SSq}])}{a(29)}$$

$$\mathrm{joint\_coeff} = \sqrt{w_{\rm VDS} \times w_{\rm DVP}}$$

$$\mathrm{variant\_branch} = |w_{\rm VDS} - w_{\rm DVP}|$$

The geometric mean $\sqrt{w_{\rm VDS} \times w_{\rm DVP}}$ is the **field coupling** between
the continuous-vacuum (VDS) and discrete-vortex (DVP) representations of the same underlying
DPM vacuum geometry.

### 5.2 Triple-Point Condition

At $[\mathrm{SSq}] = 0.57$, simultaneous convergence (within 1\%) holds:

$$\left|\frac{g_C - g_R}{g_R}\right| < 0.01, \quad \left|\frac{g_R - g_{\rm BSFG}}{g_{\rm BSFG}}\right| < 0.01, \quad \left|\frac{g_C - g_{\rm BSFG}}{g_{\rm BSFG}}\right| < 0.01$$

where $g_C$ = Compressed-MUGE, $g_R$ = Resonance-MUGE, $g_{\rm BSFG}$ = buoyancy gravity.
This triple convergence is what determines $[\mathrm{SSq}] = 0.57$ uniquely
(see PAPER\_1154 for derivation).

---

## 6. BH26$\times$BSH Cross-Resonance

The BH26$\times$BSH resonance evaluates the Buoyancy Harmonic Series (PAPER\_429) at each
spectral frequency bin of the BH26 ladder:

$$f_k = \frac{f_{\rm RR}}{\lambda_k} \quad [\mathrm{Hz}]$$

$$\mathrm{resonance}_k = \mathrm{BSH}(\omega_k) \cdot \cos(\pi t_n)$$

$$\mathrm{energy\_density}_k = \mathrm{resonance}_k \times \rho_{\rm SCm} \quad [\mathrm{J/m}^3]$$

This cross-resonance provides the **frequency-domain bridge** between the discrete BH26
spectral ladder and the continuous SCm vacuum density field.

---

## 7. Test Suite: T61--T80

| Test | Branch | Assertion |
|------|--------|-----------|
| T61 | VDS | $\mathrm{vds\_prime} > 0$ |
| T62 | VDS | $\mathrm{vds\_density} > 0$ |
| T63 | VDS | $\mathrm{vds\_k\_weighted} > 0$ |
| T64 | DVP | $\mathrm{zeta\_sum} > 0$ |
| T65 | DVP | $n_{\rm primes} > 0$ |
| T66 | DVP | $a(29) > 0$ |
| T67 | BH26 | $\Sigma_{10} = 1760$ |
| T68 | BH26 | $\mathrm{casimir\_energy} > 0$ |
| T69 | BH26 | $\mathrm{degeneracy} = 26$ |
| T70 | BH26 | $\mathrm{vds\_coupling} > 0$ |
| T71 | VDS$\times$DVP | $\mathrm{joint\_coeff} > 0$ |
| T72 | VDS$\times$DVP | $\mathrm{variant\_branch} \geq 0$ |
| T73 | VDS$\times$DVP | $w_{\rm VDS} > 0$ |
| T74 | VDS$\times$DVP | $w_{\rm DVP} > 0$ |
| T75 | BSH Res. | $\mathrm{energy\_density} \geq 0$ |
| T76 | BSH Res. | $\mathrm{freq\_k} > 0$ |
| T77 | BSH Res. | $|\mathrm{resonance}| \geq 0$ |
| T78 | Coupled | $g_{\rm BSFG} = a_{\rm DPM} \times \mathrm{joint\_coeff} \times \mathrm{vds\_k\_weighted}$ |
| T79 | Convergence | $|\mathrm{vds\_prime} - 1| < 0.1$ |
| T80 | Integration | All branch results non-zero simultaneously |

**Result: 80/80 tests PASS** (commits 6e69e2f5, dad3ae0d)

---

## 8. CP4 Calculator Classes

Session 202 produced five CP4 Phase-4 calculator classes (CP4 \#644--\#648):

| \# | Class | PAPER |
|----|-------|-------|
| 644 | `VDSBranchCalculator` | PAPER\_1151 |
| 645 | `DVPBranchCalculator` | PAPER\_1151 |
| 646 | `BH26BranchCalculator` | PAPER\_1151 |
| 647 | `VDSDVPCoupledCalculator` | PAPER\_1151 |
| 648 | `BH26BSHResonanceCalculator` | PAPER\_1151 |

---

## 9. Conclusions

Phase H202 completes the branch-level derivation of all three UQFF number systems (VDS, DVP, BH26)
and their coupled-field geometry. The geometric-mean joint coupling coefficient
$\sqrt{w_{\rm VDS} \times w_{\rm DVP}}$ is the fundamental UQFF calibration quantity for
the BSFG gravity path. All results are verified at $[\mathrm{SSq}] = 0.57$.

## References

1. Murphy, D.T. (2026). *QCalcGeom.py v2.1.0 Phase H202.* Star-Magic repository, commit 6e69e2f5.
2. Murphy, D.T. (2026). *Star-Magic.txt Chapter 4: Three Canonical Number Systems.*
3. Murphy, D.T. (2025). *PAPER\_429: BSH Buoyancy Harmonic Series.*
4. Murphy, D.T. (2025). *PAPER\_533: VDS-DVP-BH Catalogue Hub.* Session 143.
5. Murphy, D.T. (2025). *PAPER\_535: VDS-DVP-BH extended.* Session 145.
6. Kaluza, T. (1921). Zum Unitätsproblem der Physik. *Sitzungsber. Preuss. Akad. Wiss.* 966--972.
7. Klein, O. (1926). Quantentheorie und fünfdimensionale Relativitätstheorie. *Z. Phys.* 37:895--906.

*Updated: 2026-05-06 (Session 234, PAPER\_1151). Compliant: CVW v2.0.0, G6 SM Anchor Gate.
Author: Daniel T. Murphy*
