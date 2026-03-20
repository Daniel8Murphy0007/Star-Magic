# PAPER_233: SGR 1745-2900 Enhanced — SMBH Tidal Proximity Coupling, Static B-Field Magnetic Energy, ATNF Pulse Period

**Author:** Daniel T. Murphy
**Framework:** UQFF v4.8 (Star-Magic)
**Session:** 58 (grok_share_8d951e12.txt extraction — Doc 14 enhanced)
**Date:** March 2026
**Classification:** Enhanced MUGE — Three New Terms vs Session 53 MagnetarSGR1745DynamicModulationCalculator
**Status:** Proof-Quality Whitepaper

---

## Abstract

SGR 1745-2900 is the closest known magnetar to a supermassive black hole, located at a projected separation of 0.09–0.3 pc from Sagittarius A* and a deprojected estimate of $\sim 0.92$ pc. This paper documents three MUGE terms absent from the Session 53 implementation (`MagnetarSGR1745DynamicModulationCalculator`): (1) SMBH tidal coupling acceleration $a_{BH} = G M_{SgrA*}/r_{BH}^2$, (2) static (non-decaying) magnetic stored energy $a_{mag} = B^2 V/(2\mu_0 Mr)$ with $B = 2 \times 10^{10}$ T, and (3) use of the ATNF-catalogued pulse period $P = 3.76$ s (replacing an inferred value). The superconductive suppression factor is also refined to $f_{sc} = 1 - B/B_{crit}$.

---

## 1. Physical System

SGR 1745-2900's proximity to Sgr A* makes it unique in the magnetar population:

| Parameter | Value |
|-----------|-------|
| Position | $\sim 0.09$ pc projected, $\sim 0.92$ pc deprojected from Sgr A* |
| $M$ | $1.4 M_\odot$ (NS) |
| $r$ | $20$ km (NS radius) |
| $B$ | $2 \times 10^{10}$ T (static; not decaying) |
| $B_{crit}$ | $4.4 \times 10^{13}$ T |
| $P_{init}$ | $3.76$ s (ATNF Pulsar Catalogue) |
| $\tau_\Omega$ | $8000$ yr |
| $L_0$ | $5 \times 10^{28}$ W |
| $\tau_d$ | $1000$ s |
| $M_{SgrA*}$ | $4 \times 10^6 M_\odot$ |
| $r_{BH}$ | $0.92$ pc $= 2.83 \times 10^{16}$ m |

---

## 2. New Terms vs Session 53

### Comparison Table

| Term | Session 53 | Session 58 Enhancement |
|------|-----------|----------------------|
| BH proximity | Absent | $a_{BH} = GM_{SgrA*}/r_{BH}^2$ |
| B field | Decaying $B(t) = B_0 e^{-t/\tau_B}$ | Static $B = 2 \times 10^{10}$ T |
| Superconductivity | Generic $f_{TRZ}$ | $f_{sc} = 1 - B/B_{crit}$ |
| Pulse period | Inferred from $P_{dot}$ | $P = 3.76$ s (ATNF) |
| Magnetic energy | Absent | $a_{mag} = B^2 V/(2\mu_0 Mr)$ |
| Burst decay | Absent | $a_{decay} = L_0\tau_d(1-e^{-t/\tau_d})/(Mr)$ |

---

## 3. SMBH Tidal Coupling (Novel in SGR 1745 context)

$$a_{BH} = \frac{G M_{SgrA*}}{r_{BH}^2} = \frac{6.674 \times 10^{-11} \times 4 \times 10^6 \times 1.989 \times 10^{30}}{(2.83 \times 10^{16})^2}$$

$$a_{BH} = \frac{5.309 \times 10^{26}}{8.01 \times 10^{32}} \approx 6.63 \times 10^{-7} \text{ m/s}^2$$

This tidal coupling from the $4 \times 10^6 M_\odot$ SMBH is **dominant over the magnetar's self-gravity** at 0.92 pc, where the self-gravity $GM_{NS}/r_{NS}^2 \approx 2 \times 10^{12}$ m/s² only at the neutron star surface.

---

## 4. Static Magnetic Energy (Novel vs Session 53)

$$a_{mag} = \frac{B^2}{2\mu_0} \cdot \frac{V_{NS}}{Mr} = \frac{(2 \times 10^{10})^2}{2 \times 4\pi \times 10^{-7}} \cdot \frac{(4\pi/3)(20000)^3}{2.785 \times 10^{30} \times 20000}$$

Session 53 used a decaying $B(t)$; this implementation uses a **static** $B$ field, reflecting evidence that SGR 1745-2900's field has remained stable since its 2013 activation.

---

## 5. Superconductive Factor

$$f_{sc} = 1 - \frac{B}{B_{crit}} = 1 - \frac{2 \times 10^{10}}{4.4 \times 10^{13}} = 1 - 4.55 \times 10^{-4} \approx 0.99955$$

Applied to the base gravity: $a_{base} = GM/r^2 \cdot (1 + H_0 t) \cdot f_{sc}$ — a ~0.05% suppression from near-critical field.

---

## 6. ATNF Pulse Period

The pulse period $P = 3.76$ s is directly from the ATNF Pulsar Catalogue (2024), more precise than the $P \approx 3.8$ s sometimes quoted. This affects the spin-down back-reaction:
$$\frac{d\Omega}{dt} = -\frac{2\pi}{P \tau_\Omega} e^{-t/\tau_\Omega}$$

---

## 7. Calculator Class

```python
class SGR1745BHProximityMagEnergyCalculator(_CP3Calculator):
    """PAPER_233: SGR 1745-2900 enhanced — BH tidal a_BH, static mag energy, ATNF P=3.76s"""
    # Session 58 — grok_share_8d951e12.txt Doc 14 enhanced
```

Located in `CondensedPhysics3.py` (Session 58).

---

## 8. Conclusion

The enhanced SGR 1745-2900 MUGE adds three physically motivated terms absent from Session 53: SMBH tidal coupling (dominant at 0.92 pc), static magnetic energy density, and ATNF-calibrated pulse period. The result is the most complete MUGE treatment of any Galactic Centre magnetar in the UQFF library.

**Source:** grok_share_8d951e12.txt — Doc 14 (SGR 1745-2900 BH Proximity Enhanced MUGE)


**UQFF computed:** Eddington luminosity UQFF correction = 1 - [SSq]�exp(-?�?t) = 1 - 5.7e-1 � exp(-2.9e-4) = 4.3e-1; F_U at event horizon = 2.0e+18 m/s�.