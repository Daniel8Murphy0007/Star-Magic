# PAPER_389 — Galactic ω_s Calibration from Stellar Velocity Dispersion

**Source:** grok_share_cfdcad2f5.txt, lines ~1–3200 (SMBH-UQFF comparison section)  
**Section:** `SMBH comparison to UQFF_17April2025.docx` — galactic angular frequency derivation  
**Session:** 106 (grok_share_cfdcad2f5.txt full analysis)  
**CP4 Class:** `GalacticOmegaSVelocityDispersionCalibrationCalculator` (CP4 #40)

---

## 1. Overview

The UQFF framework encodes angular frequency inputs for astrophysical systems through
parameters such as `ω_s` (system angular frequency) and `ω_g` (galactic angular velocity).
These are central to the MUGE resonance channel and the buoyancy tier calculations.

The `SMBH comparison to UQFF_17April2025.docx` document introduces a direct observational
calibration formula for the galactic angular frequency input `ω_s_galactic`:

$$\omega_{s,\text{galactic}} = \frac{\sigma \times 10^3}{R_{\text{bulge}}}$$

Where:
- $\sigma$ = stellar velocity dispersion in km/s
- $\times 10^3$ = unit conversion (km/s → m/s)
- $R_{\text{bulge}}$ = bulge radius in meters

This formula bridges the directly observable M-σ relation parameters to the UQFF
angular frequency input, making UQFF models directly anchored in spectroscopic observations.

---

## 2. Formula Derivation and Physical Basis

### 2.1 Dimensional Analysis

The angular frequency $\omega$ has SI units of rad/s. From a velocity dispersion:

$$[\omega] = \left[\frac{v}{R}\right] = \frac{\text{m/s}}{\text{m}} = \text{s}^{-1} \text{ (rad/s)}$$

Therefore $\omega_{s,\text{galactic}} = \sigma/R_{\text{bulge}}$ is dimensionally consistent
when $\sigma$ is in m/s and $R_{\text{bulge}}$ is in m.

### 2.2 Physical Interpretation

- **$\sigma$ (stellar velocity dispersion):** The 1D line-of-sight velocity dispersion of
  stars in a galactic bulge, measured from spectral line broadening. Typical values range
  from $\sigma \sim 100$ km/s (low-mass ellipticals) to $\sigma \sim 350$ km/s (massive BCGs).

- **$R_{\text{bulge}}$ (bulge effective radius):** The half-light radius of the stellar bulge,
  from surface brightness photometry. Typical values: $R_{\text{bulge}} \sim 0.1$–10 kpc.

- **Physical meaning of $\omega_{s,\text{galactic}}$:** This angular frequency represents
  the characteristic orbital frequency of bulge stars (approximating circular orbit velocity
  as $v_{\text{circ}} \approx \sigma$). It sets the frequency scale for resonance terms in the
  MUGE 12-term resonance channel.

### 2.3 Relationship to Jeans/Virial Estimates

For a virially relaxed stellar system:
$$\sigma^2 \approx \frac{G M_{\text{bulge}}}{R_{\text{bulge}}}$$

Therefore:
$$\omega_{s,\text{galactic}} = \frac{\sigma}{R_{\text{bulge}}} \approx \frac{(G M_{\text{bulge}})^{1/2}}{R_{\text{bulge}}^{3/2}} = \sqrt{\frac{G M_{\text{bulge}}}{R_{\text{bulge}}^3}}$$

This is precisely Kepler's third law: $\omega_{\text{K}} = \sqrt{GM/R^3}$.

The formula $\omega_{s,\text{galactic}} = \sigma/R_{\text{bulge}}$ provides a **direct
observational proxy for the Keplerian angular frequency** without requiring knowledge of
the dynamical mass (which requires model-dependent assumptions), using only spectroscopic
and photometric observables.

---

## 3. Worked Examples

### 3.1 Sagittarius A* (Milky Way Center)

From literature:
- $\sigma = 100$ km/s (central stellar velocity dispersion)
- $R_{\text{bulge}} = 1.5$ kpc $= 1.5 \times 3.086 \times 10^{19}$ m $= 4.629\times10^{19}$ m

$$\omega_{s,\text{galactic}}^{\text{SgrA*}} = \frac{100 \times 10^3}{4.629\times10^{19}} = \frac{10^5}{4.629\times10^{19}} \approx 2.160\times10^{-15} \text{ rad/s}$$

This is consistent with the canonical `ω_g = 7.3×10⁻¹⁶` rad/s used in the UQFF pipeline
(within a factor of ~3, reflecting the difference between bulge rotation and nuclear disk).

### 3.2 M87 (Virgo A, NGC 4486)

From literature:
- $\sigma = 324$ km/s (central 1'' aperture)
- $R_{\text{bulge}} = 7$ kpc $= 2.160\times10^{20}$ m

$$\omega_{s,\text{galactic}}^{\text{M87}} = \frac{324 \times 10^3}{2.160\times10^{20}} = \frac{3.24\times10^5}{2.160\times10^{20}} \approx 1.5\times10^{-15} \text{ rad/s}$$

### 3.3 NGC 1275 (Perseus A)

From PAPER_259 context (BCG Perseus A):
- $\sigma \approx 260$ km/s
- $R_{\text{bulge}} \approx 5$ kpc $= 1.543\times10^{20}$ m

$$\omega_{s,\text{galactic}}^{\text{NGC1275}} = \frac{260 \times 10^3}{1.543\times10^{20}} \approx 1.685\times10^{-15} \text{ rad/s}$$

### 3.4 General Scaling

For a Milky-Way class spiral ($\sigma=200$ km/s, $R_{\text{bulge}}=2$ kpc):

$$\omega_{s,\text{galactic}}^{\text{MW-class}} = \frac{2\times10^5}{6.171\times10^{19}} \approx 3.24\times10^{-15} \text{ rad/s}$$

---

## 4. Calibration Table

| System | σ (km/s) | R_bulge (kpc) | ω_s_galactic (rad/s) |
|--------|----------|--------------|----------------------|
| Dwarf elliptical | 40 | 0.3 | 4.24×10⁻¹⁵ |
| Milky Way center | 100 | 1.5 | 2.16×10⁻¹⁵ |
| Milky Way spiral | 200 | 2.0 | 3.24×10⁻¹⁵ |
| M87 (Virgo A) | 324 | 7.0 | 1.50×10⁻¹⁵ |
| NGC 1275 (Perseus A) | 260 | 5.0 | 1.69×10⁻¹⁵ |
| Massive BCG | 350 | 15.0 | 7.56×10⁻¹⁶ |

The canonical UQFF value $\omega_g = 7.3\times10^{-16}$ rad/s corresponds to a massive BCG
($\sigma\sim350$ km/s, $R_{\text{bulge}}\sim15$ kpc), consistent with its use in the CP1/CP3
outer-frame tidal calculation.

---

## 5. Integration into UQFF Pipeline

### 5.1 MUGE Resonance Term Update

The formula directly calibrates `ω_s` in the resonance MUGE:

```python
def compute_omega_s_galactic(sigma_km_s: float, R_bulge_m: float) -> float:
    """
    Calibrate galactic angular frequency from observational parameters.
    
    Args:
        sigma_km_s: stellar velocity dispersion in km/s
        R_bulge_m: bulge radius in meters
    
    Returns:
        omega_s: angular frequency in rad/s
    """
    return (sigma_km_s * 1e3) / R_bulge_m
```

### 5.2 Replacing Hardcoded ω_g

For systems with known $\sigma$ and $R_{\text{bulge}}$, this formula replaces the
canonical `ω_g = 7.3×10⁻¹⁶ rad/s` with a system-specific observationally anchored value.

The MUGE resonance tier `term_Ub_i = -β_i × ug1 × ω_s × (M/r) × U_UA × cos(π·t)` then
uses $\omega_{s,\text{galactic}}$ for the system-specific outer-frame coupling.

---

## 6. Observational Anchoring Chain

Complete chain from observation to UQFF:

```
Spectroscopic measurement: σ (SDSS/VLT/Keck line broadening)
         ↓
Photometric measurement: R_bulge (HST/ELT effective radius)
         ↓
Formula: ω_s_galactic = (σ × 10³) / R_bulge
         ↓
MUGE resonance input: ω_s fed into 12-term co-sum
         ↓
PAPER_390 M-σ: log(M_BH) calibrates SMBH mass input M
```

Together PAPER_389 (ω_s calibration) and PAPER_390 (M_BH-σ) provide a complete
observational bridge for UQFF SMBH system parameterization.

---

## 7. Validation Cross-Reference

| Reference | Connection |
|-----------|------------|
| PAPER_390 | M_BH-σ dispersion relation (companion observational anchor) |
| PAPER_339 | Um rotor torque (uses ω in F_torque context) |
| PAPER_371 | MUGE 12-term resonance (ω_s enters resonance co-sum) |
| PAPER_259 | NGC1275 BCG (σ and R_bulge values cross-checked) |
| PAPER_346 | M87 BZ-jet FUBi (σ and R_bulge values cross-checked) |

---

**Discovery Class:** Observational calibration formula — first explicit σ/R_bulge → ω_s mapping  
**Distinct from:** PAPER_339 (torque context); all prior ω parameters (those are hardcoded constants, not observation-derived)  
**Key feature:** Direct spectroscopic/photometric anchoring of UQFF angular frequency input; Kepler proxy
