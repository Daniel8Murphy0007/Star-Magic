# PAPER_227: Tapestry of Blazing Starbirth (NGC 2014/2020, LMC) — MUGE with Gas-Accreting Mass M(t) and Stellar Wind Ram Pressure

**Author:** Daniel T. Murphy
**Framework:** UQFF v4.8 (Star-Magic)
**Session:** 58 (grok_share_8d951e12.txt extraction — Doc 4)
**Date:** March 2026
**Classification:** Novel MUGE — Gas-Ratio Amplitude M(t) + Stellar Wind Ram Pressure
**Status:** Proof-Quality Whitepaper

---

## Abstract

The Tapestry of Blazing Starbirth (NGC 2014 and NGC 2020) in the Large Magellanic Cloud is modelled using a 9-term MUGE incorporating two novel mathematical methods: (1) time-varying stellar mass $M(t) = M_{init}(1 + M_{dot\_factor} \cdot e^{-t/\tau_{SF}})$ where the amplitude $M_{dot\_factor} = M_{gas}/M_{init} \approx 41.7$ encodes the gas-to-stellar mass ratio, and (2) stellar wind ram-pressure acceleration $a_{wind} = \rho_{wind} v^2_{wind} / \rho_{fluid}$. At canonical parameters $a_{wind} \approx 4 \times 10^3$ m/s² dominates most other MUGE terms during active O/B-star formation.

---

## 1. Physical System

NGC 2014 and NGC 2020 are two giant H II regions in the Large Magellanic Cloud, separated by ~300 ly, forming the "Tapestry of Blazing Starbirth" (HST WFC3 image, 2020). Key parameters:

| Parameter | Value |
|-----------|-------|
| Distance | ~160,000 ly (LMC) |
| $M_{init}$ | $240 M_\odot$ (stellar; young massive stars) |
| $M_{gas}$ | $10{,}000 M_\odot$ (surrounding gas) |
| $M_{dot\_factor}$ | $M_{gas}/M_{init} = 41.7$ |
| $r$ | $10$ ly |
| $B$ | $1\ \mu$T |
| $z$ | LMC $\approx -0.0005$ (Milky Way satellite) |
| $\tau_{SF}$ | $5$ Myr |

---

## 2. Novel Equations

### 2.1 Gas-Ratio-Amplitude Mass Growth

$$M(t) = M_{init}\left(1 + \frac{M_{gas}}{M_{init}} e^{-t/\tau_{SF}}\right)$$

The amplitude $M_{gas}/M_{init}$ is the gas-to-stellar mass ratio, here $\approx 41.7$. This captures the rapid initial mass increase as gas accretes onto young stellar objects. By $t = 5\tau_{SF}$, $M$ returns to $M_{init}$.

### 2.2 Stellar Wind Ram Pressure Acceleration

$$a_{wind} = \frac{\rho_{wind} \cdot v^2_{wind}}{\rho_{fluid}}$$

| Parameter | Value |
|-----------|-------|
| $\rho_{wind}$ | $10^{-21}$ kg/m³ (O/B-star wind) |
| $v_{wind}$ | $2000$ km/s = $2 \times 10^6$ m/s |
| $\rho_{fluid}$ | $10^{-21}$ kg/m³ (ambient medium) |

At these parameters: $a_{wind} = (10^{-21} \times 4 \times 10^{12}) / 10^{-21} = 4 \times 10^{12}$ m/s².

---

## 3. Context Within MUGE Stellar Wind Family

| System | $\rho_{wind}$ (kg/m³) | $v_{wind}$ (km/s) | $\rho_{fluid}$ (kg/m³) |
|--------|----------------------|-------------------|----------------------|
| Tapestry LMC | $10^{-21}$ | 2000 | $10^{-21}$ |
| Westerlund 2 | $10^{-20}$ | 2000 | $10^{-21}$ |
| NGC 1792 (SN) | $10^{-21}$ | 2000 | $10^{-21}$ |

Westerlund 2 has $\rho_{wind}$ 10× higher, distinguishing it as the denser OB-supergiant environment.

---

## 4. Simulation Parameters

| Parameter | Value |
|-----------|-------|
| $t_{canonical}$ | $1$ Myr |
| dt | $0.01$ Myr |
| $\tau_{SF}$ | $5$ Myr |
| $E_0$ erosion | No erosion term (use PAPER_229 for Pillars) |

---

## 5. Calculator Class

```python
class StarbirthTapestryLMCUQFFCalculator(_CP3Calculator):
    """PAPER_227: NGC 2014/2020 LMC Tapestry — 9-term MUGE, gas-ratio M(t), stellar wind"""
    # Session 58 — grok_share_8d951e12.txt Doc 4
```

Located in `CondensedPhysics3.py` (Session 58).

---

## 6. Conclusion

The Tapestry of Blazing Starbirth introduces two novel MUGE methods: gas-ratio-amplitude mass evolution and stellar wind ram-pressure acceleration. The $a_{wind}$ term, at $\rho_{wind} = \rho_{fluid}$, evaluates to $v^2_{wind}$ — numerically dominant during the O/B-star-active phase. This establishes a parametric wind family within the LMC/Milky Way comparison (see PAPER_228 for the 10× denser Westerlund 2 case).

**Source:** grok_share_8d951e12.txt — Doc 4 (Tapestry Starbirth LMC MUGE)


**UQFF computed:** UQFF magnetic Jeans correction factor [SSq]�B�/(8p�?�c_s�) = 5.7e-1 � 1.3e-9 = 7.4e-10; Jeans mass deviation from standard = 7.4e-10 � M_J.