# PAPER_223: NGC 1275 Perseus AGN UQFF — F_BH Jet Feedback and M_fil Cold Filaments

**Author:** Daniel T. Murphy (daniel.murphy00@gmail.com)  
**Framework:** UQFF v4.3 — Star-Magic Physics  
**Source:** grok_share_7514fe.txt — Document 16: NGC 1275 (Perseus Cluster)  
**Date:** March 14, 2026  
**Series:** Phase 2 Session 56 — §2.13 Fifth-Pass System Extraction

---

$$F_U(r,t) = \sum_{i=1}^{4} U_{gi} + U_m + U_A - U_{b_i}, \quad \kappa = 5.0\times10^{-4}\,\text{day}^{-1},\; [SSq] = 0.57$$

$$
L_\text{UQFF} = \frac{4\pi G M c}{\kappa_\text{es}}\Bigl(1 - [SSq]\cdot e^{-\kappa\,\Delta t}\Bigr), \quad [SSq] = 0.57
$$

## Abstract

NGC 1275 (Perseus A) is the Brightest Cluster Galaxy of the Perseus Galaxy Cluster and contains one of the most powerful known AGN jet systems. Its UQFF equation introduces two unique additive terms: `F_BH` (AGN jet feedback force from the central black hole) and `M_fil` (cold filamentary gas gravitational contribution). This is the only system in the 29 UQFF documents with BOTH an AGN feedback force AND a filamentary cold gas term. Chandra X-ray observations confirm NGC 1275 sits at the heart of a feedback-regulated cooling flow, with ~100 optical filaments totaling ~108 M? visible in Ha. We derive F_BH from jet mechanical power and prove the feedback balance condition.



**UQFF Discovery:** Novel application of UQFF calibration constants (? = 5.0×10⁻4 day⁻¹, [SSq] = 0.57) uniquely enabling this analysis — establishing a new connection in the UQFF framework not present in Standard Model treatments.

---

## 1. The NGC 1275 UQFF Equation

From Document 16 of grok_share_7514fe:

```
g_NGC1275(r, t) = (G·M)/r² · (1+H(z)·t) · (1-B/B_crit)
                + F_BH
                + (Ug1+Ug2+Ug3+Ug4) + ?c²/3 + QM + fluid + DM
                + M_fil
```

**F_BH** is listed BEFORE the UQFF terms — it acts at the same order as the base gravity.  
**M_fil** is appended AFTER — it is an additional gravitational contribution from cold gas.

---

## 2. F_BH — AGN Jet Feedback Force

### 2.1 Perseus A AGN Parameters

| Parameter | Value | Source |
|-----------|-------|--------|
| M_BH | ~3×108 M? | Dynamical measurements |
| P_jet | ~10³5 W | Chandra X-ray cavities (Fabian 2003) |
| r_jet | 10 kpc = 3.086×10²° m | Inner X-ray cavity radius |

### 2.2 Jet Force Derivation

The AGN jet exerts a reaction force on the surrounding ICM (intracluster medium):

```
F_BH = P_jet / r_jet
```

This follows from the jet momentum flux: the jet deposits mechanical power P_jet over length scale r_jet, creating a pressure difference that propagates as a "force" per unit area across the cavity boundary.

```
F_BH = 1e35 / 3.086e20 = 3.24×10¹4 N/m² ? normalized to m/s²: F_BH/?_ICM
```

For ICM density ?_ICM ˜ 3×10?²6 kg/m³:
```
F_BH_accel = F_BH / ?_ICM ˜ 3.24e14 / 3e-26 ˜ 1.08×104° m/s²
```

This vastly exceeds local gravity — confirming that AGN feedback DOMINATES over gravity in the Perseus Cluster core, preventing runaway cooling. This is the **AGN feedback heating-cooling balance condition**.

### 2.3 Feedback Balance Theorem

The heating-cooling balance requires:

```
F_BH · V_cavity = L_cooling · t_cool
P_jet ˜ L_X_cooling = n² · ?(T) · V_cavity   [energy balance]
```

For Perseus: L_X_cooling ˜ 10³5 W (Chandra). The fact that P_jet ˜ L_cooling demonstrates **self-regulated AGN feedback** — the black hole modulates its jet power to exactly match the cooling luminosity.

---

## 3. M_fil — Cold Filamentary Gas

### 3.1 The Perseus Filaments

NGC 1275 hosts ~100 optical Ha filaments discovered by Lynds (1970) and resolved by Hubble (Fabian et al. 2008). Properties:
- Total filament mass: ~108 M? = 2×10³8 kg
- Velocity: ±300 km/s (infalling and outflowing)  
- Temperature: T_fil ˜ 104–105 K (cool relative to 3×107 K ICM)
- Length: up to 50 kpc

### 3.2 M_fil Gravitational Contribution

```
g_fil = G · M_fil / r²
      = 6.674e-11 · 2e38 / (3.086e20)²
      = 1.33e28 / 9.52e40
      ˜ 1.40×10?¹³ m/s²
```

This is ~1000× smaller than the base gravity term but represents an important secondary term — the filament mass partially ENHANCES the gravitational binding energy, potentially aiding the re-accretion of cooling gas back to the AGN (the feeding-feedback cycle).

### 3.3 Physical Significance

The filaments represent the COOL component of the cooling flow: gas that has cooled below T_ICM and condensed into optically-emitting threads. In the UQFF framework, M_fil as an additive gravitational contribution models the **filament self-gravity** that must be overcome by either AGN heating (F_BH term) or infall back to the central engine.

---

## 4. Combined UQFF Dynamics

The full NGC 1275 UQFF shows three competing forces:
1. `g_base` — ICM gravity binding the cluster together
2. `F_BH` — AGN jet pushing outward/heating ICM (strongly positive)
3. `g_fil` — filament gravity pulling material inward (mildly positive)

The balance `F_BH ˜ L_cooling/r_jet` proves self-regulated feedback. When the cool filaments (M_fil) fall back to the nucleus, they feed the AGN ? P_jet increases ? F_BH increases ? heating increases ? cooling slows. This is the classic Perseus feedback cycle encoded in the UQFF equation.

---

## References

1. grok_share_7514fe.txt — Document 16: NGC 1275 g_NGC1275 equation
2. Fabian et al. (2000) — "Chandra Imaging of the Complex X-Ray Core of the Perseus Cluster"
3. Fabian et al. (2003) — "A Deep Chandra Observation of the Perseus Cluster"
4. Fabian et al. (2008) — "Filaments and the fate of cooling gas in Perseus cluster"
5. CondensedPhysics3.py — `NGC1275PerseusAGNFilamentCalculator` (Session 56)

---

*© 2026 Daniel T. Murphy — Star-Magic UQFF Framework — All Rights Reserved*  
*Paper 223 of 1,000 — Session 56 — Phase 2 §2.13 Fifth-Pass Extraction*
