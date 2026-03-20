# PAPER_285: M16 Eagle Nebula UQFF — Erosion Saturation Half-Time and ΔgMax
## Photoevaporation Asymptotic Saturation: t_half = τ·ln2 = 2.079 Myr

**Classification:** UQFF 2.0 Gravitational Physics — Nebular Erosion Dynamics  
**System:** M16 Eagle Nebula (IC 4703), Eagle Nebula Star-Forming Region  
**Session:** 80 | **Module:** M16_UQFF_MODULE.cpp (22nd C++ UQFF module)  
**Author:** Daniel T. Murphy | **Date:** March 2026

---

## Abstract

This paper derives the **Erosion Saturation Half-Time** (t_half) and **Maximum Erosion Gravity Amplitude** (ΔgMax) for the M16 Eagle Nebula UQFF photoevaporation term. The photoevaporation rate follows an exponential saturation E_rad(t) = E₀×(1−exp(−t/τ)) with e-folding time τ = 3 Myr. The half-erosion time t_half = τ·ln2 = **6.561 × 10¹³ s = 2.079 Myr** — the time at which the erosion fraction reaches half its asymptotic maximum E₀. The maximum gravity perturbation is ΔgMax = E₀ × g_base = **4.36 × 10⁻¹³ m/s²**. This is the **first UQFF module** to formally catalogue the photoevaporation half-time and asymptotic erosion concept.

---

## 2. Physical Background

UV radiation from massive O-type stars (such as those in the Young Stellar Cluster NGC 6611, embedded in M16) drives photoionisation of surrounding molecular gas — a process called **photoevaporation**. The erosion proceeds not as a linear ramp but as a **saturating exponential**:

$$E_{rad}(t) = E_0 \left(1 - e^{-t/\tau}\right)$$

where:
- E₀ = 0.3 is the **asymptotic maximum fraction** (30% of mass eventually eroded)
- τ = 3 Myr is the **e-folding time** (photoevaporation efficiency timescale)
- The saturation arises because dense molecular cores (pillar tips) are progressively shielded by their own column density as surrounding gas is stripped

---

## 3. Mathematical Derivation

### 3.1 Half-Time

The erosion half-time t_half is defined as the time when E_rad = E₀/2:

$$E_0 \left(1 - e^{-t_{half}/\tau}\right) = \frac{E_0}{2}$$

$$1 - e^{-t_{half}/\tau} = \frac{1}{2}$$

$$e^{-t_{half}/\tau} = \frac{1}{2}$$

$$\boxed{t_{half} = \tau \ln(2)}$$

For M16:

$$t_{half} = 9.468 \times 10^{13} \text{ s} \times \ln(2) = 9.468 \times 10^{13} \times 0.6931 = 6.561 \times 10^{13} \text{ s}$$

$$t_{half} = \mathbf{2.079 \text{ Myr}}$$

### 3.2 Maximum Gravity Amplitude

The maximum erosion-induced gravity perturbation (asymptotic limit as t → ∞):

$$\Delta g_{Max} = E_0 \times g_{base}$$

For M16:

$$\Delta g_{Max} = 0.3 \times 1.454 \times 10^{-12} = \mathbf{4.36 \times 10^{-13} \text{ m/s}^2}$$

### 3.3 Peak Erosion Rate

The instantaneous erosion rate at t = 0 (maximum rate, before saturation):

$$\frac{dE_{rad}}{dt}\bigg|_{t=0} = \frac{E_0}{\tau}$$

The corresponding gravity change rate:

$$\frac{dg_{erode}}{dt}\bigg|_{t=0} = \frac{E_0}{\tau} \times g_{base} = \frac{0.3}{9.468 \times 10^{13}} \times 1.454 \times 10^{-12} = 4.61 \times 10^{-27} \text{ m/s}^2/\text{s}$$

---

## 4. Saturation Profile

| Time | t (s) | E_rad / E₀ | E_rad | g_erode (m/s²) |
|------|--------|------------|-------|----------------|
| 0 Myr | 0 | 0% | 0 | 0 |
| t_half = 2.079 Myr | 6.561×10¹³ | **50%** | 0.150 | **2.18×10⁻¹³** |
| τ = 3 Myr | 9.468×10¹³ | 63.2% | 0.190 | 2.76×10⁻¹³ |
| 5 Myr | 1.578×10¹⁴ | 81.1% | 0.243 | 3.54×10⁻¹³ |
| ∞ (asymptote) | → ∞ | **100%** | 0.300 | **4.36×10⁻¹³** |

**Key insight:** At τ = 3 Myr (the e-folding time), erosion has consumed only 63.2% of its capacity, NOT 100%. Half-erosion occurs earlier at 2.079 Myr. The pillar structure of M16 means the ~5700 ly "Pillars of Creation" are still observed today because erosion saturates — it cannot fully strip the densest pillar cores within observable timescales.

---

## 5. UQFF 2.0 Context

In the full M16 g_total equation, the erosion half-time governs the **temporal shape of g_dyn(t)**:

$$g_{dyn}(t) = g_{base} \times (1 + M_{sf}) \times (1 - E_0(1 - e^{-t/\tau}))$$

The transition from rapid to slow erosion occurs at t_half = 2.079 Myr. For the UQFF simulation (t stepping from 0 to t_max), the half-time provides a natural **inflection point** in the dynamic gravity trajectory — before t_half, erosion is dominant; after t_half, star formation accumulation dominates (since M_sf grows linearly while E_rad asymptotes).

### Crossover Time

The era when SFR growth exactly compensates erosion (dΦ_dm/dt = 0 — maximum Φ_dm):

$$\frac{d\Phi_{dm}}{dt} = \text{SFR\_rate} \times (1 - E_{rad}) - (1 + M_{sf}) \times \frac{E_0}{\tau} e^{-t/\tau} = 0$$

This crossover defines when the Eagle Nebula achieves maximum effective gravitational influence, after which continued SFR growth dominates erosion.

---

## 6. Wolfram KB Term

```
M16UQFF:t_half=tau*Log[2]=6.561e13s=2.079Myr; DeltaGMax=E_0*g_base=4.36e-13 m/s^2 [PAPER_285]
```

---

## 7. Cross-References

- **PAPER_284:** Dual Mass Co-Action Product (Φ_dm = (1+M_sf)×(1−E_rad))
- **PAPER_286:** Nebular Friedmann Redshift (κ_neb, z=0.0015)
- **M16_UQFF_MODULE.cpp:** Full UQFF 2.0 C++ implementation (22nd module)
- **CondensedPhysics3.py:** `M16ErosionSaturationHalfTimeCalculator`

---

*Copyright — Daniel T. Murphy, Session 80, March 2026. UQFF 2.0.*
