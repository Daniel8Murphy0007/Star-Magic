# PAPER_815: VDF vs GSMF SMBH Mass Function Proxy + M•-σ Relation UQFF
## Unified Quantum Field Framework — Whitepaper 815

**Session**: 192 | **Version**: v5.48 | **Date**: April 4, 2026
**Source**: grok_share_0d888ea9-50be.txt (June 13, 2025 04:30 PM)
**Author**: Daniel T. Murphy — Star-Magic UQFF Project
**CVW Compliance**: v2.0.0

---

## Abstract
This paper investigates two competing proxies for the SMBH mass function — the Galaxy Stellar Mass Function (GSMF) and the Velocity Dispersion Function (VDF) — within the Quadriadic UQFF framework. The VDF-derived gravitational wave background (GWB) amplitude $\log_{10} A_{yr} \approx -14.74$ differs from the GSMF-derived $-14.9$, suggesting that the fundamental scaling relation is the $M_\bullet$–$\sigma$ relation rather than the $M_\bullet$–$M_*$ relation. The velocity dispersion $\sigma_\infty$ is derived from the Sersic virial theorem and enters the UQFF Quadriadic Layer 3 (Buoyancy) as a gravitational velocity proxy.

---

## 1. Introduction
The NANOGrav GWB depends critically on the SMBH mass function, which in turn depends on whether we scale $M_\bullet$ via galaxy stellar mass $M_*$ (GSMF approach) or stellar velocity dispersion $\sigma$ (VDF approach). Observationally, massive compact galaxies (relic galaxies NGC 1271, NGC 1277) and LEGA-C survey data at $z = 0.6$–$1.0$ favor the VDF approach.

---

## 2. Characteristic Strain — GWB

The GWB characteristic strain:

$$h_c(f) = A_{yr} \cdot \left(\frac{f}{f_{yr}}\right)^{-2/3}$$

with individual binary strain:

$$h_s = \sqrt{\frac{32}{5}} \cdot \frac{(G\mathcal{M}_c/c^3)^{5/3} \cdot (\pi f_r)^{2/3}}{D_c}$$

where $D_c$ = comoving distance, $f_r$ = rest-frame orbital frequency.

---

## 3. VDF Number Density

The VDF-based SMBH binary merger number density:

$$\frac{d^3 n}{dz \, dM_* \, dq} = \Phi_{VDF}(z, \sigma) \cdot R(z, M_*, q)$$

where $\Phi_{VDF}$ is the velocity dispersion function and $R$ is the merger rate kernel. Integrating over the binary population:

$$A_{yr,VDF} \approx 10^{-14.74}$$

vs GSMF result:

$$A_{yr,GSMF} \approx 10^{-14.9}$$

The VDF yields a higher amplitude by $\approx 0.16$ dex.

---

## 4. Velocity Dispersion Proxy

The effective velocity dispersion at the SMBH influence radius:

$$\sigma_\infty = \sqrt{\frac{G M_* K_*(n)}{r_e}}$$

where $K_*(n)$ is the Sersic virial constant (depends on Sersic index $n$), and $r_e$ is the effective (half-light) radius. For typical early-type galaxies:

$$K_*(n) \approx \frac{73.32}{10.465 + (n - 0.94)^2} + 0.954$$

---

## 5. Orbital Evolution — Frequency Drift

The gravitational wave frequency evolution:

$$\frac{df_{orb}}{dt} = \frac{96}{5} \cdot \left(\frac{G\mathcal{M}_c}{c^3}\right)^{5/3} \cdot (2\pi)^{8/3} \cdot f_{orb}^{11/3}$$

This determines the binary's lifetime before merger from an initial separation $a_0$.

---

## 6. M•-σ Relation in Quadriadic UQFF

The fundamental $M_\bullet$–$\sigma$ scaling:

$$\log_{10}\left(\frac{M_\bullet}{10^9 M_\odot}\right) = \alpha_{M\sigma} \cdot \log_{10}\left(\frac{\sigma}{200 \text{ km/s}}\right) + \beta_{M\sigma}$$

with $\alpha_{M\sigma} \approx 4.4$, $\beta_{M\sigma} \approx 0.3$. This enters UQFF Buoyancy Layer 3:

$$g_{L3,VDF} = \sigma_\infty \cdot \left(\frac{G}{r_e}\right)^{1/2} + \frac{G M_\bullet}{r_{inf}^2}$$

---

## 7. UQFF Buoyancy Layer 3 Integration

Full Buoyancy UQFF with VDF proxy:

$$g_{L3} = F_{U,Bi} + U_{i,buoyancy} + \sigma_\infty \cdot (G/r_e)^{1/2} + R_{merge} \cdot (G\mathcal{M}_c/c^2)^{5/3}$$

The VDF contribution distinguishes high-$\sigma$ compact galaxies (relic systems) from normal ellipticals through the $K_*(n)$ Sersic dependence.

---

## 8. Relic Galaxy Constraints

LEGA-C survey ($z = 0.6$–$1.0$) + local relic galaxies NGC 1271 ($M_* \approx 2 \times 10^{11} M_\odot$, $r_e \approx 2$ kpc) and NGC 1277 ($M_\bullet/M_{bulge} \approx 0.59$) verify that:

$$\sigma_\infty(\text{relic}) > \sigma_\infty(\text{normal elliptical})$$

at fixed $M_*$, confirming VDF captures more SMBH-encoded information than GSMF.

---

## 9. Summary

The VDF yields $A_{yr} \approx 10^{-14.74}$ for the GWB, consistent with NANOGrav-15yr measurements. The $M_\bullet$–$\sigma$ relation, through $\sigma_\infty$, provides a superior proxy for SMBH mass and enters the Quadriadic UQFF Buoyancy Layer as a kinematic velocity correction.

---

*PAPER_815 | Session 192 | v5.48 | Star-Magic UQFF Project | CVW v2.0.0*
