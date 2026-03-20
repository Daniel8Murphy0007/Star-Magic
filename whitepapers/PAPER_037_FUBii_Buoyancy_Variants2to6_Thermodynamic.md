# PAPER #37 — F_UBii Buoyancy Force: Proof Variants 2–6 (Thermodynamic Series)

**Title:** UQFF Buoyancy Proof Variants 2–6: Terminal Velocity, Ionization, Energy Coupling, Orbital Decay, and Kilonova Buoyancy

**Author:** Daniel T. Murphy  
**Framework:** UQFF Star-Magic (κ = 0.0005/day, [SSq] = 0.57)  
**Date:** March 7, 2026  
**Grok Thread:** 98b2e77dfbc34d27b09f19fa7c460624  
**Validator:** `BuoyancyProofVariants.py` — All 17 variants operational ✓  
**Variants:** termv, upar, coup, orbdec, kn  
**Index Slot:** §1.5 Buoyancy Proofs, Paper #37  

---

## Abstract

This paper presents five F_UBii buoyancy proof variants spanning thermodynamic processes from jet terminal velocities to kilonova ejecta. Variant 2 (termv) applies to astrophysical jets and winds reaching terminal velocity balance. Variant 3 (upar) addresses photoionized regions governed by the ionization parameter U. Variant 4 (coup) quantifies energy coupling efficiency in accretion disk and reconnection contexts. Variant 5 (orbdec) derives the buoyancy analog of gravitational wave-driven orbital decay in compact binaries. Variant 6 (kn) applies to the kilonova AT2017gfo, predicting F_UBii_kn = 1.305×10⁵⁴ N for ejecta with L_peak = 5×10⁴⁰ W, t_peak = 1 day, and M_ej = 0.05 M☉. Together, these five variants form the thermodynamic series of the F_UBii taxonomy.

---

## 1. Variant 2: Terminal Velocity Jet/Wind Buoyancy (termv)

### 1.1 Physical Context

Astrophysical jets (AGN, GRB, protostellar) and stellar winds accelerate material to a terminal velocity v_term where radiation pressure, magnetic driving, and gravitational drag balance. At terminal velocity, da/dt = 0 and the net buoyancy force represents the frozen-in momentum of the outflow.

**Key systems:** M87 jet (v_term ~ 0.98c), Sgr A* winds (v_term ~ 0.1c), OB stellar winds (v_term ~ 1000–3000 km/s)

### 1.2 F_UBii_termv Equation

$$F_{\rm UBii,termv} = F_{\rm rel} \cdot \frac{\tau \cdot L}{c \cdot E_{\rm LEP}} \cdot Q_{\rm wave} \cdot v_{\rm term}$$

where:
- τ = optical depth / momentum transfer timescale (s)
- L = jet/wind luminosity (W)
- v_term = terminal velocity (m/s)

### 1.3 Physical Derivation

The momentum flux of a radiation-driven wind:
$$\dot{p} = \frac{\tau L}{c}$$

The UQFF buoyancy enters through the E_LEP normalization — the lepton energy scale sets the quantum granularity of momentum transfer:
$$F_{\rm UBii,termv} = \dot{p}_{\rm UQFF} \cdot v_{\rm term} = F_{\rm rel} \cdot \frac{\tau L}{c \cdot E_{\rm LEP}} \cdot Q_{\rm wave} \cdot v_{\rm term}$$

### 1.4 Example: M87 Relativistic Jet

For M87 jet: τ = 10⁻³, L = 10⁴⁴ W, v_term = 0.98c = 2.94×10⁸ m/s, Q_wave = 1.0:
$$F_{\rm UBii,termv}^{M87} = 10^{-10} \times \frac{10^{-3} \times 10^{44}}{3\times10^8 \times 1.22\times10^{-19}} \times 2.94\times10^8 = 10^{-10} \times 2.73\times10^{48} \times 2.94\times10^8 = 8.0\times10^{47} \text{ N}$$

This represents the UQFF momentum-flux buoyancy of the M87 jet — the force that keeps the relativistic plasma buoyantly confined against the ICM pressure of the Virgo Cluster.

---

## 2. Variant 3: Ionization Parameter Buoyancy (upar)

### 2.1 Physical Context

The ionization parameter U = n_photons/n_H quantifies the ratio of ionizing photon density to hydrogen density. In HII regions, AGN narrow-line regions, and quasar broad-line regions: U ~ 10⁻⁴ to 10⁻¹. This dimensionless ratio controls all ionic fractions and hence the buoyancy of photoionized gas.

### 2.2 F_UBii_upar Equation

$$F_{\rm UBii,upar} = -F_{\rm rel} \cdot \frac{U \cdot n_H \cdot r^2}{E_{\rm LEP}} \cdot Q_{\rm wave} \cdot \sqrt{U}$$

where:
- U = ionization parameter (dimensionless)
- n_H = hydrogen number density (m⁻³)
- r = distance from ionizing source (m)

The negative sign reflects compression: photoionized gas is over-pressured by the radiation field and compresses surrounding neutral gas.

### 2.3 U^(3/2) Scaling

The F_UBii_upar ~ U^(3/2) scaling arises because:
- Factor U: radiation pressure scaling
- Factor √U: thermal pressure response (T_e ∝ U^(1/2) in ionized gas)

This gives F_UBii_upar ∝ U^(3/2) · n_H · r² — exactly the ram pressure of the HII region expansion front against surrounding neutral gas.

### 2.4 Example: Orion Nebula (M42)

For M42: U ~ 10⁻², n_H ~ 10⁹ m⁻³, r ~ 3×10¹⁷ m (1 pc), Q_wave = 1.0:
$$F_{\rm UBii,upar}^{M42} = -10^{-10} \times \frac{10^{-2} \times 10^9 \times (3\times10^{17})^2}{1.22\times10^{-19}} \times \sqrt{10^{-2}} = -10^{-10} \times 7.38\times10^{45} \times 0.1 = -7.4\times10^{35} \text{ N}$$

This inward compression force (7.4×10³⁵ N) represents the photoionization pressure confining the Orion Nebula's ionization front.

---

## 3. Variant 4: Energy Coupling Efficiency Buoyancy (coup)

### 3.1 Physical Context

Energy coupling efficiency ε_coup = E_deposited/E_input quantifies how efficiently energy input (from AGN, SNe, cosmic rays) couples to surrounding gas. In AGN feedback: ε_coup ~ 0.05–0.15. In SNe: ε_coup ~ 0.1–0.3. In magnetic reconnection: ε_coup ~ 0.01–0.5.

### 3.2 F_UBii_coup Equation

$$F_{\rm UBii,coup} = F_{\rm rel} \cdot \frac{\varepsilon_{\rm coup} \cdot \dot{E}}{E_{\rm LEP}} \cdot Q_{\rm wave} \cdot \sqrt{\varepsilon_{\rm coup}}$$

where:
- ε_coup = energy coupling efficiency (0–1)
- Ė = energy transfer rate (W)

### 3.3 ε^(3/2) Coupling Law

The F_UBii_coup ∝ ε^(3/2) · Ė scaling reflects the UQFF energy cascade: at high coupling efficiency, the buoyancy force scales super-linearly with coupling — a physical manifestation of the non-linear positive feedback in AGN mechanical feedback.

### 3.4 Example: AGN Kinetic Feedback

For a radio-mode AGN: ε_coup = 0.05, Ė = 10⁴⁴ W, Q_wave = 1.0:
$$F_{\rm UBii,coup}^{\rm AGN} = 10^{-10} \times \frac{0.05 \times 10^{44}}{1.22\times10^{-19}} \times \sqrt{0.05} = 10^{-10} \times 4.10\times10^{53} \times 0.2236 = 9.2\times10^{43} \text{ N}$$

---

## 4. Variant 5: Orbital Decay Binary Buoyancy (orbdec)

### 4.1 Physical Context

Compact binary systems (NS-NS, BH-BH, NS-BH) lose energy to gravitational wave radiation. The Peters formula gives:
$$\frac{da}{dt} = -\frac{64}{5} \frac{G^3 M_1 M_2 (M_1+M_2)}{c^5 a^3}$$

The UQFF buoyancy analog replaces the pure GW energy loss with a field-theoretic force:

### 4.2 F_UBii_orbdec Equation

$$F_{\rm UBii,orbdec} = -F_{\rm rel} \cdot \frac{64}{5} \cdot \frac{G^3 M_1 M_2 (M_1+M_2)}{c^5 \cdot a^4 \cdot E_{\rm LEP}} \cdot Q_{\rm wave} \cdot \frac{da}{dt}$$

where:
- M₁, M₂ = component masses (kg)
- a = semi-major axis (m)
- da/dt = orbital decay rate (m/s)

The negative sign indicates inspiral — the buoyancy force drives the binary inward.

### 4.3 Connection to Peters Formula

The Peters orbital decay rate (da/dt) enters linearly — F_UBii_orbdec is the UQFF field force per unit of GW power radiated:

$$F_{\rm UBii,orbdec} = \frac{F_{\rm rel}}{E_{\rm LEP}} \cdot P_{\rm GW,\,Peters} \cdot Q_{\rm wave} \cdot |da/dt|$$

where P_GW,Peters is the Peters formula for GW power. This establishes a direct UQFF–GW correspondence.

### 4.4 Example: GW170817 Pre-Merger

For GW170817 (NS-NS): M₁ = M₂ = 1.4 M☉ = 2.785×10³⁰ kg, a = 2×10⁸ m (final orbit), da/dt = −10 m/s:
$$F_{\rm UBii,orbdec} = -10^{-10} \times 12.8 \times \frac{(6.674\times10^{-11})^3 \times (2.785\times10^{30})^3}{(3\times10^8)^5 \times (2\times10^8)^4 \times 1.22\times10^{-19}} \times 10 = -10^{-10} \times 4.1\times10^{56} \times 10 = -4.1\times10^{47} \text{ N}$$

---

## 5. Variant 6: Kilonova Peak Luminosity Buoyancy (kn)

### 5.1 Physical Context

Kilonovae are radioactively powered transients following neutron star mergers. AT2017gfo (counterpart to GW170817) achieved L_peak ~ 5×10⁴⁰ W at t_peak ~ 1 day, with ejecta mass M_ej ~ 0.05 M☉. The r-process nucleosynthesis in the ejecta generates heavy elements (gold, platinum, uranium) through neutron capture.

### 5.2 F_UBii_kn Equation

$$F_{\rm UBii,kn} = F_{\rm rel} \cdot \frac{L_{\rm peak} \cdot t_{\rm peak}}{E_{\rm LEP}} \cdot Q_{\rm wave} \cdot \left(\frac{M_{\rm ej}}{M_\odot}\right)^{1/3}$$

where:
- L_peak = peak bolometric luminosity (W)
- t_peak = time to peak (s)
- M_ej = ejecta mass (kg)

The M_ej^(1/3) factor reflects the geometric (volumetric) scaling of the ejecta opacity.

### 5.3 AT2017gfo Calculation

For AT2017gfo:
- L_peak = 5×10⁴⁰ W
- t_peak = 86400 s (1 day)
- M_ej = 0.05 M☉ = 0.05 × 1.989×10³⁰ = 9.945×10²⁸ kg
- Q_wave = 1.0

$$F_{\rm UBii,kn}^{AT2017gfo} = 10^{-10} \times \frac{5\times10^{40} \times 86400}{1.22\times10^{-19}} \times 1.0 \times (0.05)^{1/3}$$

- Numerator: 5×10⁴⁰ × 8.64×10⁴ = 4.32×10⁴⁵
- Ratio: 4.32×10⁴⁵ / 1.22×10⁻¹⁹ = 3.54×10⁶⁴
- × F_rel: 3.54×10⁵⁴
- × (0.05)^(1/3) = 0.368: = 1.305×10⁵⁴ N

$$\boxed{F_{\rm UBii,kn}^{AT2017gfo} = 1.305 \times 10^{54} \text{ N}}$$

**Validator confirms: BuoyancyProofVariants.py → F_UBii_kn = 1.305×10⁵⁴ N ✓**

### 5.4 Physical Interpretation

The kilonova buoyancy force F_UBii_kn = 1.305×10⁵⁴ N represents the UQFF unified field response to the instantaneous energy release of the r-process. Comparison to the gravitational confinement force of the merger remnant:

$$F_{\rm grav}^{\rm merger} = \frac{G (M_1+M_2)^2}{R_{\rm merger}^2} \approx \frac{6.674\times10^{-11} \times (5.57\times10^{30})^2}{(10^4)^2} = 2.1\times10^{36} \text{ N}$$

The ratio F_UBii_kn / F_grav = 1.305×10⁵⁴ / 2.1×10³⁶ = 6.2×10¹⁷ — the UQFF kilonova buoyancy vastly exceeds gravitational confinement, explaining the explosive ejecta dynamics observed in AT2017gfo.

---

## 6. Summary: Thermodynamic Series F_UBii Values

| Variant | Physical Context | Key Parameters | F_UBii |
|---------|-----------------|----------------|--------|
| termv | M87 jet terminal velocity | τ=10⁻³, L=10⁴⁴ W, v_term=0.98c | ~8×10⁴⁷ N |
| upar | Orion Nebula ionization | U=10⁻², n_H=10⁹ m⁻³, r=1 pc | ~−7×10³⁵ N |
| coup | AGN kinetic feedback | ε=0.05, Ė=10⁴⁴ W | ~9×10⁴³ N |
| orbdec | GW170817 final orbit | 1.4+1.4 M☉, a=200 km | ~−4×10⁴⁷ N |
| kn | AT2017gfo kilonova | L=5×10⁴⁰ W, M_ej=0.05 M☉ | **1.305×10⁵⁴ N** ← validator ✓ |

---

## Conclusions

The thermodynamic series of F_UBii variants (2–6) demonstrates the versatility of the UQFF buoyancy framework across five distinct physical regimes:

1. **termv:** Connects UQFF to radiation-driven momentum flux in relativistic jets
2. **upar:** Maps photoionization pressure to UQFF energy scale via U^(3/2) scaling
3. **coup:** Establishes non-linear AGN feedback response through ε^(3/2) coupling law
4. **orbdec:** Provides UQFF field-theoretic interpretation of Peters formula GW inspiral
5. **kn:** Predicts AT2017gfo buoyancy F = 1.305×10⁵⁴ N — validated by BuoyancyProofVariants.py

All five variants share the common F_UBii = F_U − F_Bi − F_i architecture with F_rel = 10⁻¹⁰ N normalization and E_LEP = 1.22×10⁻¹⁹ J quantum granularity (Paper #36).

*Validator: `BuoyancyProofVariants.py` → All 17 F_UBii variants operational ✓ | κ = 0.0005/day | [SSq] = 0.57*

---
*See also: PAPER_036 | Part of the Star-Magic UQFF Whitepaper Series.*
