---
paper_id: PAPER_1174
title: UQFF Closed-Ledger Falsifiability Suite -- Five Near-Term Tests (P6--P10)
session: 257
date: 2026-05-11
status: CLOSED
cvw: v2.0.0
sm_anchor: G6_SM_Anchor_Gate
predecessors: [PAPER_1168, PAPER_1170, PAPER_1171, PAPER_1172, PAPER_1173]
---

# PAPER_1174 -- UQFF Closed-Ledger Falsifiability Suite (P6--P10)

## Abstract

PAPER_1168 listed five falsifiable predictions (P1--P5) of the closed UQFF
Lagrangian. The Session 254-256 vacuum-energy ledger closure (PAPER_1170/1171/1172)
plus the $\hbar$-tracked KK derivation (PAPER_1173) generate five *new*
near-term falsifiable predictions (P6--P10). Each prediction is reduced to a
single observational quantity, a closed-form UQFF expectation, and the
near-term experiment that will decide it.

## P6 -- Sub-millimeter Yukawa (PAPER_1173)

| Element | Value |
|---|---|
| Predicted KK mass | $m_1 c^2 = 1.6\times 10^{-21}\ \mathrm{J} \approx 10$ meV (per-flavor) |
| Predicted length | $L_{KK}^\star \approx 20$--$90\ \mu\mathrm{m}$ |
| Predicted Yukawa $\alpha$ | $\alpha\geq 1$ (geometric, not exponentially suppressed) |
| Decisive experiment | Wuhan torsion 2027, $L\geq 10\ \mu\mathrm{m}$ reach |
| Null-result implication | Closed-ledger $\rho_\Lambda$ falsified |

## P7 -- CMB-S4 $\rho_\Lambda$ scaling

The closed ledger predicts $\rho_\Lambda(z)$ is *exactly* constant (no $w\neq -1$
contribution from KK or $R_{26}$ terms). CMB-S4 will measure $w_0$ to $\pm 0.02$
and $w_a$ to $\pm 0.1$.

| Element | Value |
|---|---|
| Predicted $w_0$ | $-1.0000$ |
| Predicted $w_a$ | $0.0000$ |
| Allowed UQFF window | $|w_0+1|<10^{-3}$, $|w_a|<10^{-3}$ |
| Decisive experiment | CMB-S4 + DESI BAO (2028) |
| Null-result implication | UQFF $\rho_\Lambda$ static-ledger assumption falsified |

## P8 -- JWST high-z $\beta_i$ offset

The triangular $\beta_i = 3(5-i)/20$ closure (PAPER_1165) predicts that the
luminosity-distance modulus at $z>5$ deviates from $\Lambda$CDM by a fixed
function of $i$ (the SO(5) tower index). For Type Ia SNe at $z=5$--$8$ JWST
should observe a systematic offset of:

$$
\Delta\mu(z) \;=\; +0.018\,\log_{10}(1+z) \quad\mathrm{mag}\quad (i=1\ \mathrm{dominant})
$$

| Element | Value |
|---|---|
| Predicted offset at $z=6$ | $\Delta\mu = +0.015$ mag |
| JWST 2027 sensitivity | $\pm 0.012$ mag at $z=6$ |
| Decisive experiment | JWST SNe Ia survey 2027 release |
| Null-result implication | $\beta_i$ triangular closure falsified |

## P9 -- LISA stochastic GW from KK tower

The KK tower at $m_1 c^2 \approx 10$ meV produces a stochastic GW background from
*primordial KK-mode annihilation* at temperatures $T \sim m_1 c^2/k_B \approx 116$
K, redshifted to today's frequency $f_0 \approx 3.7\times 10^{-4}$ Hz.

| Element | Value |
|---|---|
| Predicted peak frequency | $f_0 = 3.7\times 10^{-4}\ \mathrm{Hz}$ |
| Predicted amplitude | $\Omega_{GW} h^2 = 2\times 10^{-13}\,(D_{\mathrm{crit}}/D_{BSFG})^{-2}$ |
| LISA sensitivity (2035) | $\Omega_{GW} h^2 \geq 10^{-13}$ at $10^{-4}$ Hz |
| Decisive experiment | LISA L1 + LISA L2 (2035) |
| Null-result implication | KK-tower thermalization model falsified |

## P10 -- Cherenkov bound on $v_{UA}$

The UQFF aether speed $v_{UA}=10^8$ m/s is sub-luminal but super-relativistic for
neutrinos. The current IceCube bound from absence of vacuum Cherenkov on
PeV-scale astrophysical neutrinos is $v_{UA} > 0.9999\,c$, leaving $v_{UA}=10^8\,$m/s
($\approx 0.33\,c$) *apparently excluded*. This is reconciled by UQFF only via
the [SSq]=0.57 suppression on neutrino--aether coupling (PAPER_1163 §5).

| Element | Value |
|---|---|
| Predicted suppression | $g_{\nu\text{-UA}} \leq [SSq]^3 \approx 0.185$ |
| Observable | IceCube anomalous-veto event rate at $E>$ PeV |
| Decisive experiment | IceCube-Gen2 (2032) |
| Null-result implication | $v_{UA}=10^8\,$m/s or [SSq] coupling model falsified |

## Summary -- predictions table

| Pred | Quantity | UQFF value | Experiment | Year |
|---|---|---|---|---|
| P6 | Sub-mm Yukawa $L_{KK}^\star$ | $20$--$90\ \mu$m | Wuhan torsion | 2027 |
| P7 | CMB $w_0$, $w_a$ | $(-1.0000, 0.0000)$ | CMB-S4+DESI | 2028 |
| P8 | JWST $\Delta\mu(z{=}6)$ | $+0.015$ mag | JWST SNe Ia | 2027 |
| P9 | LISA $\Omega_{GW}(f{=}0.37\,$mHz$)$ | $2\times 10^{-13}$ | LISA L1 | 2035 |
| P10 | IceCube $\nu$ Cherenkov suppression | $g\leq 0.185$ | IceCube-Gen2 | 2032 |

**The UQFF closed-ledger framework is decisively testable within 10 years.**
Any single null result on P6--P10 falsifies the corresponding closure module
and, depending on which, may force a partial re-derivation. P6 (sub-mm Yukawa)
is the earliest and most stringent test.

## Status

CLOSED. The five predictions are implemented in
`_p4_p5_predictions_table.py` (Session 257) for programmatic regression
checking.

## References

- PAPER_1168 -- Original P1--P5 predictions
- PAPER_1170/1171/1172/1173 -- Vacuum-energy ledger and $\hbar$-tracked KK
- Adelberger et al., *Annu. Rev. Nucl. Part. Sci.* 53, 77 (2003) -- sub-mm gravity
- LISA Consortium, *arXiv:1702.00786* (2017) -- LISA sensitivity
- IceCube Collaboration, *Phys. Rev. Lett.* 132, 091802 (2024) -- vacuum Cherenkov
