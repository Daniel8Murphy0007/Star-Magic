# PAPER_816: EHT, ngEHT, and BHEX New SMBH Population Photon Ring UQFF
## Unified Quantum Field Framework — Whitepaper 816

**Session**: 192 | **Version**: v5.48 | **Date**: April 4, 2026
**Source**: grok_share_0d888ea9-50be.txt (June 13, 2025 04:45 PM); "Accessing_a_New_Population_of_Supermassive_Black_H.pdf" (ApJ May 2025)
**Author**: Daniel T. Murphy — Star-Magic UQFF Project
**CVW Compliance**: v2.0.0

---

## Abstract
This paper integrates the May 2025 ApJ study "Accessing a New Population of Supermassive Black Holes" into the Quadriadic UQFF framework. Twelve new SMBH targets ($\log M_\bullet/M_\odot = 8.7$–$9.7$) are resolvable at the photon-ring scale by the ngEHT (Next-Generation Event Horizon Telescope) and the proposed BHEX (Black Hole Explorer) space-based extension. The photon ring angular diameter, Eddington accretion ratio $f_{Edd}$, electron-to-ion temperature ratio $T_e/T_i$, and Blandford-Znajek (BZ) jet mechanism parameters all enter the UQFF Quadriadic layers.

---

## 1. Introduction
Following the historic imaging of M87* (2019) and Sgr A* (2022), the EHT and its successors aim to image a broader SMBH population. The photon ring angle:

$$\theta_{ring} = 2\sqrt{27} \cdot \frac{G M_\bullet}{c^2 D}$$

provides a direct mass-distance measurement. The May 2025 ApJ study identifies 12 new SMBH targets in the photon-ring accessible regime with ngEHT (resolution ~15 Gλ).

---

## 2. New SMBH Targets — 12 Systems

| Source | $\log M_\bullet/M_\odot$ | $D$ (Mpc) | $\theta_{ring}$ (µas) | $f_{Edd}$ |
|--------|-------------------------|-----------|----------------------|----------|
| IC 1459 | 9.45 | 29 | 8.9 | ~10⁻³ |
| NGC 4594 (Sombrero) | 8.82 | 9.8 | — | ~10⁻⁴ |
| NGC 4261 | 8.72 | 31 | — | ~10⁻⁴ |
| NGC 315 | 9.0 | 70 | — | 1.7×10⁻⁴ |
| NGC 1218 | 8.7 | 50 | — | ~10⁻³ |
| 7 additional | 8.7–9.7 | varied | — | ~10⁻⁴–10⁻³ |

All are Low-Luminosity AGN (LLAGN) in RIAF (Radiatively Inefficient Accretion Flow) regime.

---

## 3. Photon Ring Angular Diameter

The critical impact parameter of the photon ring for Schwarzschild BH:

$$\theta_{ring} = 2\sqrt{27} \cdot \frac{G M_\bullet}{c^2 D}$$

For IC 1459 ($M_\bullet = 2.8 \times 10^9 M_\odot$, $D = 29$ Mpc):

$$\theta_{ring} = 2\sqrt{27} \cdot \frac{(6.67\times10^{-11})(5.6\times10^{39})}{(9\times10^{16})^2 \cdot (8.96\times10^{23})} \approx 8.9 \; \mu\text{as}$$

This enters UQFF Layer 1 as:

$$g_{L1,ring} = g_{UQFF}(r,t) + \left(\frac{G M_\bullet}{c^2 D}\right) \cdot \theta_{ring}$$

---

## 4. Eddington Ratio and RIAF

The Eddington accretion ratio:

$$f_{Edd} = \frac{\dot{M}}{\dot{M}_{Edd}} = \frac{L_{bol}}{\eta c^2} \cdot \frac{\kappa_{es} c}{4\pi G M_\bullet}$$

For $f_{Edd} \ll 0.01$, the accretion is in the RIAF (hot, geometrically thick) regime. This enters Layer 2:

$$g_{L2,Edd} = f_{Edd} \cdot (T_e/T_i) \cdot \frac{G M_\bullet}{r^2}$$

---

## 5. Electron-to-Ion Temperature Ratio

In the RIAF, the electron and ion temperatures differ strongly:

$$\frac{T_e}{T_i} = R_{high} \cdot \frac{\beta^2}{1+\beta^2} + R_{low} \cdot \frac{1}{1+\beta^2}$$

where $\beta = P_{gas}/P_{mag}$, $R_{high} \approx 40$ (weakly magnetized), $R_{low} \approx 1$. For typical $\beta = 1$:

$$T_e/T_i = 40 \cdot 0.5 + 1 \cdot 0.5 = 20.5$$

---

## 6. Blandford-Znajek (BZ) Mechanism

The jet power from a spinning BH in an external magnetic field:

$$P_{BZ} = \frac{\kappa}{4\pi c} \cdot \Phi_{BH}^2 \cdot \Omega_H^2$$

where $\Phi_{BH}$ = magnetic flux threading the horizon, $\Omega_H = a_* c/(2 r_H)$ = horizon angular velocity, $\kappa \approx 0.044$ (parabolic field). This enters Layer 4:

$$g_{L4,BZ} = P_{BZ} \cdot \frac{1}{r^3} \cdot a_*$$

---

## 7. Full UQFF Integration — EHT Sources

For NGC 315 ($M_\bullet = 10^9 M_\odot$, $f_{Edd} = 1.7\times10^{-4}$, $a_* = 0.9$):

$$g_{L1} = g_{UQFF} + \theta_{ring} \cdot (G M_\bullet/c^2 D)$$
$$g_{L2} = f_{Edd} \cdot (T_e/T_i) = 1.7\times10^{-4} \times 20.5 \approx 3.5\times10^{-3}$$
$$g_{L3} = F_{U,Bi} + U_{i,buoyancy}$$
$$g_{L4} = P_{BZ}/r^3 \cdot a_*$$

---

## 8. ngEHT and BHEX Resolution

- **ngEHT**: ~15 Gλ at 1.3 mm (EHT ×5 baselines), resolves $\theta \geq$  a few µas
- **BHEX**: Space-based, planned 350+ Gλ, resolves 0.1 µas class structures
- **Photon ring substructure**: BHEX can resolve the secondary ring $n=1$ for M87*

These observational capabilities constrain $a_*$ and $M_\bullet$ to 5%–10% precision, directly calibrating the UQFF Layer 4 BZ term.

---

## 9. Summary

Twelve new SMBH targets (log M/M_sun = 8.7–9.7) are accessible to ngEHT and BHEX via photon-ring imaging. Their photon ring diameter, Eddington ratio, T_e/T_i temperature ratio, and BZ jet power all integrate directly into the Quadriadic UQFF Layers 1–4, extending the UQFF to a new population of LLAGN engines.

---

*PAPER_816 | Session 192 | v5.48 | Star-Magic UQFF Project | CVW v2.0.0*
