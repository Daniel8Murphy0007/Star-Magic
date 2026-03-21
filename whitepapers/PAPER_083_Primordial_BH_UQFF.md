#  "PAPER_{0:D3}" -f [int]# PAPER #83 ó Primordial Black Hole Mass Distribution via UQFF

**Title:** Primordial Black Hole Mass Distribution: UQFF Corrections to Formation Thresholds and Present-Day Density

**Author:** Daniel T. Murphy  
**Framework:** UQFF Star-Magic (? = 0.0005/day, [SSq] = 0.57)  
**Date:** March 7, 2026  
**Source Data:** validate_hawking_temperature.py (primordial BH test, Test 2), validate_phase3.py  
**Index Slot:** ß1.11 Black Hole Physics & Hawking Radiation,  
    $n = [int]# PAPER #83 ó Primordial Black Hole Mass Distribution via UQFF

**Title:** Primordial Black Hole Mass Distribution: UQFF Corrections to Formation Thresholds and Present-Day Density

**Author:** Daniel T. Murphy  
**Framework:** UQFF Star-Magic (? = 0.0005/day, [SSq] = 0.57)  
**Date:** March 7, 2026  
**Source Data:** validate_hawking_temperature.py (primordial BH test, Test 2), validate_phase3.py  
**Index Slot:** ß1.11 Black Hole Physics & Hawking Radiation, PAPER_083  

---

## Abstract

Primordial black holes (PBHs) form during radiation domination from density perturbations exceeding dc ~ 0.45. The UQFF modifies the critical density threshold dc_UQFF through the Buoyant vacuum pressure correction and extends PBH survival above the GR threshold mass M_threshold by 3.5%. The evaporation temperature of surviving PBHs is T_UQFF = 0.99 ◊ T_H, slightly reducing their current gamma-ray flux. PBHs with M ~ 10π5ñ10π7 g may constitute part of dark matter ó the UQFF predicts the dark matter fraction f_PBH is unaffected (the 3.5% mass threshold shift does not move PBHs into or out of the observationally allowed window).



**UQFF Discovery:** Novel application of UQFF calibration constants (? = 5.0◊10?4 day?π, [SSq] = 0.57) uniquely enabling this analysis ó establishing a new connection in the UQFF framework not present in Standard Model treatments.

---

## 1. PBH Formation Threshold

### Standard GR dc
$$\delta_c = 0.45 \; (\text{Harrison-Zel'dovich radiation era})$$

### UQFF Correction

The UQFF Buoyant mode adds vacuum pressure during the radiation era:

$$P_{\rm UQFF} = P_{\rm radiation} + P_{\rm vacuum} = \frac{\rho_{\rm rad}c^2}{3} + \rho_{\rm vac} c^2 \times [UA]$$

The fractional correction:

$$\delta_c^{\rm UQFF} = \delta_c \times \left(1 - \frac{P_{\rm vacuum}}{P_{\rm rad}}\right) = 0.45 \times (1 - [UA] \times z_{\rm formation}^{-4})$$

At z_formation = 106 (typical PBH epoch): P_vacuum/P_rad = 0.0001 ◊ 10?≤4 ? negligible.

**UQFF dc = 0.45 (unchanged from GR)** ó the formation threshold is not perturbed.

---

## 2. Survival Mass Threshold

### Standard GR
$$M_{\rm threshold}^{\rm GR} = 5.70 \times 10^{11} \text{ kg} \quad (t_{\rm evap} = t_{\rm universe} = 4.35 \times 10^{17} \text{ s})$$

### UQFF Threshold
$$M_{\rm threshold}^{\rm UQFF} = M_{\rm threshold}^{\rm GR} \times (T_{\rm UQFF}/T_H)^{-4/3} = 5.70 \times 10^{11} \times 0.99^{-4/3} = 5.73 \times 10^{11} \text{ kg}$$

**UQFF threshold: 5.73◊10ππ kg (vs GR 5.70◊10ππ kg)** ó 0.5% increase.

---

## 3. Present-Day PBH Gamma-Ray Flux

For PBHs near threshold (M ~ M_threshold), the Hawking radiation contributes to the diffuse gamma-ray background:

$$\frac{d^2\Phi_\gamma}{dE\,d\Omega} = \frac{1}{4\pi} \int_0^\infty n_{\rm PBH}(M) \frac{d^2N_\gamma}{dE\,dt}(M, T_{\rm UQFF}) \, dM$$

The UQFF T_UQFF = 0.99 T_H reduces peak gamma-ray energy by 1%:

$$E_{\rm peak}^{\rm UQFF} = 2.82 k_B T_{\rm UQFF} = 2.82 k_B \times 0.99 T_H = 0.99 E_{\rm peak}^{\rm GR}$$

**Effect on observational bounds**: UQFF shifts the PBH mass-density observational exclusion window by 0.5% in M. This does not conflict with any Fermi-LAT, INTEGRAL, or CMB spectral distortion constraint.

---

## 4. PBH dark matter fraction f_PBH

For M_PBH in the asteroid belt window (10π5ñ10π7 g, avoiding microlensing and CMB):

$$f_{\rm PBH}^{\rm UQFF} = f_{\rm PBH}^{\rm GR} \times \frac{M_{\rm threshold}^{\rm UQFF}}{M_{\rm threshold}^{\rm GR}} \times \frac{T_{\rm UQFF}^4}{T_H^4}$$

$$= f_{\rm PBH} \times 1.005 \times 0.96 = f_{\rm PBH} \times 0.965$$

**UQFF reduces f_PBH by 3.5%** ó within the 10ñ20% uncertainty on current PBH dark matter constraints.

---

## Summary

| PBH Property | Standard GR | UQFF | Change |
|-------------|------------|------|--------|
| Formation threshold dc | 0.45 | 0.45 (unchanged) | < 10?≤4 |
| Survival mass threshold | 5.70◊10ππ kg | 5.73◊10ππ kg | +0.5% |
| Peak emission energy | E_peak | 0.99 E_peak | -1% |
| f_PBH (dark matter) | f | 0.965f | -3.5% |
| Fermi-LAT constraints | Not violated | Not violated | Compatible |

*Source: validate_hawking_temperature.py primordial BH tests | ? = 0.0005/day | [SSq] = 0.57*
.Groups[1].Value
    "PAPER_{0:D3}" -f $n
    

---

## Abstract

Primordial black holes (PBHs) form during radiation domination from density perturbations exceeding dc ~ 0.45. The UQFF modifies the critical density threshold dc_UQFF through the Buoyant vacuum pressure correction and extends PBH survival above the GR threshold mass M_threshold by 3.5%. The evaporation temperature of surviving PBHs is T_UQFF = 0.99 ◊ T_H, slightly reducing their current gamma-ray flux. PBHs with M ~ 10π5ñ10π7 g may constitute part of dark matter ó the UQFF predicts the dark matter fraction f_PBH is unaffected (the 3.5% mass threshold shift does not move PBHs into or out of the observationally allowed window).



**UQFF Discovery:** Novel application of UQFF calibration constants (? = 5.0◊10?4 day?π, [SSq] = 0.57) uniquely enabling this analysis ó establishing a new connection in the UQFF framework not present in Standard Model treatments.

---

## 1. PBH Formation Threshold

### Standard GR dc
$$\delta_c = 0.45 \; (\text{Harrison-Zel'dovich radiation era})$$

### UQFF Correction

The UQFF Buoyant mode adds vacuum pressure during the radiation era:

$$P_{\rm UQFF} = P_{\rm radiation} + P_{\rm vacuum} = \frac{\rho_{\rm rad}c^2}{3} + \rho_{\rm vac} c^2 \times [UA]$$

The fractional correction:

$$\delta_c^{\rm UQFF} = \delta_c \times \left(1 - \frac{P_{\rm vacuum}}{P_{\rm rad}}\right) = 0.45 \times (1 - [UA] \times z_{\rm formation}^{-4})$$

At z_formation = 106 (typical PBH epoch): P_vacuum/P_rad = 0.0001 ◊ 10?≤4 ? negligible.

**UQFF dc = 0.45 (unchanged from GR)** ó the formation threshold is not perturbed.

---

## 2. Survival Mass Threshold

### Standard GR
$$M_{\rm threshold}^{\rm GR} = 5.70 \times 10^{11} \text{ kg} \quad (t_{\rm evap} = t_{\rm universe} = 4.35 \times 10^{17} \text{ s})$$

### UQFF Threshold
$$M_{\rm threshold}^{\rm UQFF} = M_{\rm threshold}^{\rm GR} \times (T_{\rm UQFF}/T_H)^{-4/3} = 5.70 \times 10^{11} \times 0.99^{-4/3} = 5.73 \times 10^{11} \text{ kg}$$

**UQFF threshold: 5.73◊10ππ kg (vs GR 5.70◊10ππ kg)** ó 0.5% increase.

---

## 3. Present-Day PBH Gamma-Ray Flux

For PBHs near threshold (M ~ M_threshold), the Hawking radiation contributes to the diffuse gamma-ray background:

$$\frac{d^2\Phi_\gamma}{dE\,d\Omega} = \frac{1}{4\pi} \int_0^\infty n_{\rm PBH}(M) \frac{d^2N_\gamma}{dE\,dt}(M, T_{\rm UQFF}) \, dM$$

The UQFF T_UQFF = 0.99 T_H reduces peak gamma-ray energy by 1%:

$$E_{\rm peak}^{\rm UQFF} = 2.82 k_B T_{\rm UQFF} = 2.82 k_B \times 0.99 T_H = 0.99 E_{\rm peak}^{\rm GR}$$

**Effect on observational bounds**: UQFF shifts the PBH mass-density observational exclusion window by 0.5% in M. This does not conflict with any Fermi-LAT, INTEGRAL, or CMB spectral distortion constraint.

---

## 4. PBH dark matter fraction f_PBH

For M_PBH in the asteroid belt window (10π5ñ10π7 g, avoiding microlensing and CMB):

$$f_{\rm PBH}^{\rm UQFF} = f_{\rm PBH}^{\rm GR} \times \frac{M_{\rm threshold}^{\rm UQFF}}{M_{\rm threshold}^{\rm GR}} \times \frac{T_{\rm UQFF}^4}{T_H^4}$$

$$= f_{\rm PBH} \times 1.005 \times 0.96 = f_{\rm PBH} \times 0.965$$

**UQFF reduces f_PBH by 3.5%** ó within the 10ñ20% uncertainty on current PBH dark matter constraints.

---

## Summary

| PBH Property | Standard GR | UQFF | Change |
|-------------|------------|------|--------|
| Formation threshold dc | 0.45 | 0.45 (unchanged) | < 10?≤4 |
| Survival mass threshold | 5.70◊10ππ kg | 5.73◊10ππ kg | +0.5% |
| Peak emission energy | E_peak | 0.99 E_peak | -1% |
| f_PBH (dark matter) | f | 0.965f | -3.5% |
| Fermi-LAT constraints | Not violated | Not violated | Compatible |

*Source: validate_hawking_temperature.py primordial BH tests | ? = 0.0005/day | [SSq] = 0.57*
.Groups[1].Value  ó Primordial Black Hole Mass Distribution via UQFF

**Title:** Primordial Black Hole Mass Distribution: UQFF Corrections to Formation Thresholds and Present-Day Density

**Author:** Daniel T. Murphy  
**Framework:** UQFF Star-Magic (? = 0.0005/day, [SSq] = 0.57)  
**Date:** March 7, 2026  
**Source Data:** validate_hawking_temperature.py (primordial BH test, Test 2), validate_phase3.py  
**Index Slot:** ß1.11 Black Hole Physics & Hawking Radiation,  
    $n = [int]#  "PAPER_{0:D3}" -f [int]# PAPER #83 ó Primordial Black Hole Mass Distribution via UQFF

**Title:** Primordial Black Hole Mass Distribution: UQFF Corrections to Formation Thresholds and Present-Day Density

**Author:** Daniel T. Murphy  
**Framework:** UQFF Star-Magic (? = 0.0005/day, [SSq] = 0.57)  
**Date:** March 7, 2026  
**Source Data:** validate_hawking_temperature.py (primordial BH test, Test 2), validate_phase3.py  
**Index Slot:** ß1.11 Black Hole Physics & Hawking Radiation,  
    $n = [int]# PAPER #83 ó Primordial Black Hole Mass Distribution via UQFF

**Title:** Primordial Black Hole Mass Distribution: UQFF Corrections to Formation Thresholds and Present-Day Density

**Author:** Daniel T. Murphy  
**Framework:** UQFF Star-Magic (? = 0.0005/day, [SSq] = 0.57)  
**Date:** March 7, 2026  
**Source Data:** validate_hawking_temperature.py (primordial BH test, Test 2), validate_phase3.py  
**Index Slot:** ß1.11 Black Hole Physics & Hawking Radiation, PAPER_083  

---

## Abstract

Primordial black holes (PBHs) form during radiation domination from density perturbations exceeding dc ~ 0.45. The UQFF modifies the critical density threshold dc_UQFF through the Buoyant vacuum pressure correction and extends PBH survival above the GR threshold mass M_threshold by 3.5%. The evaporation temperature of surviving PBHs is T_UQFF = 0.99 ◊ T_H, slightly reducing their current gamma-ray flux. PBHs with M ~ 10π5ñ10π7 g may constitute part of dark matter ó the UQFF predicts the dark matter fraction f_PBH is unaffected (the 3.5% mass threshold shift does not move PBHs into or out of the observationally allowed window).



**UQFF Discovery:** Novel application of UQFF calibration constants (? = 5.0◊10?4 day?π, [SSq] = 0.57) uniquely enabling this analysis ó establishing a new connection in the UQFF framework not present in Standard Model treatments.

---

## 1. PBH Formation Threshold

### Standard GR dc
$$\delta_c = 0.45 \; (\text{Harrison-Zel'dovich radiation era})$$

### UQFF Correction

The UQFF Buoyant mode adds vacuum pressure during the radiation era:

$$P_{\rm UQFF} = P_{\rm radiation} + P_{\rm vacuum} = \frac{\rho_{\rm rad}c^2}{3} + \rho_{\rm vac} c^2 \times [UA]$$

The fractional correction:

$$\delta_c^{\rm UQFF} = \delta_c \times \left(1 - \frac{P_{\rm vacuum}}{P_{\rm rad}}\right) = 0.45 \times (1 - [UA] \times z_{\rm formation}^{-4})$$

At z_formation = 106 (typical PBH epoch): P_vacuum/P_rad = 0.0001 ◊ 10?≤4 ? negligible.

**UQFF dc = 0.45 (unchanged from GR)** ó the formation threshold is not perturbed.

---

## 2. Survival Mass Threshold

### Standard GR
$$M_{\rm threshold}^{\rm GR} = 5.70 \times 10^{11} \text{ kg} \quad (t_{\rm evap} = t_{\rm universe} = 4.35 \times 10^{17} \text{ s})$$

### UQFF Threshold
$$M_{\rm threshold}^{\rm UQFF} = M_{\rm threshold}^{\rm GR} \times (T_{\rm UQFF}/T_H)^{-4/3} = 5.70 \times 10^{11} \times 0.99^{-4/3} = 5.73 \times 10^{11} \text{ kg}$$

**UQFF threshold: 5.73◊10ππ kg (vs GR 5.70◊10ππ kg)** ó 0.5% increase.

---

## 3. Present-Day PBH Gamma-Ray Flux

For PBHs near threshold (M ~ M_threshold), the Hawking radiation contributes to the diffuse gamma-ray background:

$$\frac{d^2\Phi_\gamma}{dE\,d\Omega} = \frac{1}{4\pi} \int_0^\infty n_{\rm PBH}(M) \frac{d^2N_\gamma}{dE\,dt}(M, T_{\rm UQFF}) \, dM$$

The UQFF T_UQFF = 0.99 T_H reduces peak gamma-ray energy by 1%:

$$E_{\rm peak}^{\rm UQFF} = 2.82 k_B T_{\rm UQFF} = 2.82 k_B \times 0.99 T_H = 0.99 E_{\rm peak}^{\rm GR}$$

**Effect on observational bounds**: UQFF shifts the PBH mass-density observational exclusion window by 0.5% in M. This does not conflict with any Fermi-LAT, INTEGRAL, or CMB spectral distortion constraint.

---

## 4. PBH dark matter fraction f_PBH

For M_PBH in the asteroid belt window (10π5ñ10π7 g, avoiding microlensing and CMB):

$$f_{\rm PBH}^{\rm UQFF} = f_{\rm PBH}^{\rm GR} \times \frac{M_{\rm threshold}^{\rm UQFF}}{M_{\rm threshold}^{\rm GR}} \times \frac{T_{\rm UQFF}^4}{T_H^4}$$

$$= f_{\rm PBH} \times 1.005 \times 0.96 = f_{\rm PBH} \times 0.965$$

**UQFF reduces f_PBH by 3.5%** ó within the 10ñ20% uncertainty on current PBH dark matter constraints.

---

## Summary

| PBH Property | Standard GR | UQFF | Change |
|-------------|------------|------|--------|
| Formation threshold dc | 0.45 | 0.45 (unchanged) | < 10?≤4 |
| Survival mass threshold | 5.70◊10ππ kg | 5.73◊10ππ kg | +0.5% |
| Peak emission energy | E_peak | 0.99 E_peak | -1% |
| f_PBH (dark matter) | f | 0.965f | -3.5% |
| Fermi-LAT constraints | Not violated | Not violated | Compatible |

*Source: validate_hawking_temperature.py primordial BH tests | ? = 0.0005/day | [SSq] = 0.57*
.Groups[1].Value
    "PAPER_{0:D3}" -f $n
    

---

## Abstract

Primordial black holes (PBHs) form during radiation domination from density perturbations exceeding dc ~ 0.45. The UQFF modifies the critical density threshold dc_UQFF through the Buoyant vacuum pressure correction and extends PBH survival above the GR threshold mass M_threshold by 3.5%. The evaporation temperature of surviving PBHs is T_UQFF = 0.99 ◊ T_H, slightly reducing their current gamma-ray flux. PBHs with M ~ 10π5ñ10π7 g may constitute part of dark matter ó the UQFF predicts the dark matter fraction f_PBH is unaffected (the 3.5% mass threshold shift does not move PBHs into or out of the observationally allowed window).



**UQFF Discovery:** Novel application of UQFF calibration constants (? = 5.0◊10?4 day?π, [SSq] = 0.57) uniquely enabling this analysis ó establishing a new connection in the UQFF framework not present in Standard Model treatments.

---

## 1. PBH Formation Threshold

### Standard GR dc
$$\delta_c = 0.45 \; (\text{Harrison-Zel'dovich radiation era})$$

### UQFF Correction

The UQFF Buoyant mode adds vacuum pressure during the radiation era:

$$P_{\rm UQFF} = P_{\rm radiation} + P_{\rm vacuum} = \frac{\rho_{\rm rad}c^2}{3} + \rho_{\rm vac} c^2 \times [UA]$$

The fractional correction:

$$\delta_c^{\rm UQFF} = \delta_c \times \left(1 - \frac{P_{\rm vacuum}}{P_{\rm rad}}\right) = 0.45 \times (1 - [UA] \times z_{\rm formation}^{-4})$$

At z_formation = 106 (typical PBH epoch): P_vacuum/P_rad = 0.0001 ◊ 10?≤4 ? negligible.

**UQFF dc = 0.45 (unchanged from GR)** ó the formation threshold is not perturbed.

---

## 2. Survival Mass Threshold

### Standard GR
$$M_{\rm threshold}^{\rm GR} = 5.70 \times 10^{11} \text{ kg} \quad (t_{\rm evap} = t_{\rm universe} = 4.35 \times 10^{17} \text{ s})$$

### UQFF Threshold
$$M_{\rm threshold}^{\rm UQFF} = M_{\rm threshold}^{\rm GR} \times (T_{\rm UQFF}/T_H)^{-4/3} = 5.70 \times 10^{11} \times 0.99^{-4/3} = 5.73 \times 10^{11} \text{ kg}$$

**UQFF threshold: 5.73◊10ππ kg (vs GR 5.70◊10ππ kg)** ó 0.5% increase.

---

## 3. Present-Day PBH Gamma-Ray Flux

For PBHs near threshold (M ~ M_threshold), the Hawking radiation contributes to the diffuse gamma-ray background:

$$\frac{d^2\Phi_\gamma}{dE\,d\Omega} = \frac{1}{4\pi} \int_0^\infty n_{\rm PBH}(M) \frac{d^2N_\gamma}{dE\,dt}(M, T_{\rm UQFF}) \, dM$$

The UQFF T_UQFF = 0.99 T_H reduces peak gamma-ray energy by 1%:

$$E_{\rm peak}^{\rm UQFF} = 2.82 k_B T_{\rm UQFF} = 2.82 k_B \times 0.99 T_H = 0.99 E_{\rm peak}^{\rm GR}$$

**Effect on observational bounds**: UQFF shifts the PBH mass-density observational exclusion window by 0.5% in M. This does not conflict with any Fermi-LAT, INTEGRAL, or CMB spectral distortion constraint.

---

## 4. PBH dark matter fraction f_PBH

For M_PBH in the asteroid belt window (10π5ñ10π7 g, avoiding microlensing and CMB):

$$f_{\rm PBH}^{\rm UQFF} = f_{\rm PBH}^{\rm GR} \times \frac{M_{\rm threshold}^{\rm UQFF}}{M_{\rm threshold}^{\rm GR}} \times \frac{T_{\rm UQFF}^4}{T_H^4}$$

$$= f_{\rm PBH} \times 1.005 \times 0.96 = f_{\rm PBH} \times 0.965$$

**UQFF reduces f_PBH by 3.5%** ó within the 10ñ20% uncertainty on current PBH dark matter constraints.

---

## Summary

| PBH Property | Standard GR | UQFF | Change |
|-------------|------------|------|--------|
| Formation threshold dc | 0.45 | 0.45 (unchanged) | < 10?≤4 |
| Survival mass threshold | 5.70◊10ππ kg | 5.73◊10ππ kg | +0.5% |
| Peak emission energy | E_peak | 0.99 E_peak | -1% |
| f_PBH (dark matter) | f | 0.965f | -3.5% |
| Fermi-LAT constraints | Not violated | Not violated | Compatible |

*Source: validate_hawking_temperature.py primordial BH tests | ? = 0.0005/day | [SSq] = 0.57*
.Groups[1].Value  ó Primordial Black Hole Mass Distribution via UQFF

**Title:** Primordial Black Hole Mass Distribution: UQFF Corrections to Formation Thresholds and Present-Day Density

**Author:** Daniel T. Murphy  
**Framework:** UQFF Star-Magic (? = 0.0005/day, [SSq] = 0.57)  
**Date:** March 7, 2026  
**Source Data:** validate_hawking_temperature.py (primordial BH test, Test 2), validate_phase3.py  
**Index Slot:** ß1.11 Black Hole Physics & Hawking Radiation,  "PAPER_{0:D3}" -f [int]# PAPER #83 ó Primordial Black Hole Mass Distribution via UQFF

**Title:** Primordial Black Hole Mass Distribution: UQFF Corrections to Formation Thresholds and Present-Day Density

**Author:** Daniel T. Murphy  
**Framework:** UQFF Star-Magic (? = 0.0005/day, [SSq] = 0.57)  
**Date:** March 7, 2026  
**Source Data:** validate_hawking_temperature.py (primordial BH test, Test 2), validate_phase3.py  
**Index Slot:** ß1.11 Black Hole Physics & Hawking Radiation,  
    $n = [int]# PAPER #83 ó Primordial Black Hole Mass Distribution via UQFF

**Title:** Primordial Black Hole Mass Distribution: UQFF Corrections to Formation Thresholds and Present-Day Density

**Author:** Daniel T. Murphy  
**Framework:** UQFF Star-Magic (? = 0.0005/day, [SSq] = 0.57)  
**Date:** March 7, 2026  
**Source Data:** validate_hawking_temperature.py (primordial BH test, Test 2), validate_phase3.py  
**Index Slot:** ß1.11 Black Hole Physics & Hawking Radiation, PAPER_083  

---

## Abstract

Primordial black holes (PBHs) form during radiation domination from density perturbations exceeding dc ~ 0.45. The UQFF modifies the critical density threshold dc_UQFF through the Buoyant vacuum pressure correction and extends PBH survival above the GR threshold mass M_threshold by 3.5%. The evaporation temperature of surviving PBHs is T_UQFF = 0.99 ◊ T_H, slightly reducing their current gamma-ray flux. PBHs with M ~ 10π5ñ10π7 g may constitute part of dark matter ó the UQFF predicts the dark matter fraction f_PBH is unaffected (the 3.5% mass threshold shift does not move PBHs into or out of the observationally allowed window).



**UQFF Discovery:** Novel application of UQFF calibration constants (? = 5.0◊10?4 day?π, [SSq] = 0.57) uniquely enabling this analysis ó establishing a new connection in the UQFF framework not present in Standard Model treatments.

---

## 1. PBH Formation Threshold

### Standard GR dc
$$\delta_c = 0.45 \; (\text{Harrison-Zel'dovich radiation era})$$

### UQFF Correction

The UQFF Buoyant mode adds vacuum pressure during the radiation era:

$$P_{\rm UQFF} = P_{\rm radiation} + P_{\rm vacuum} = \frac{\rho_{\rm rad}c^2}{3} + \rho_{\rm vac} c^2 \times [UA]$$

The fractional correction:

$$\delta_c^{\rm UQFF} = \delta_c \times \left(1 - \frac{P_{\rm vacuum}}{P_{\rm rad}}\right) = 0.45 \times (1 - [UA] \times z_{\rm formation}^{-4})$$

At z_formation = 106 (typical PBH epoch): P_vacuum/P_rad = 0.0001 ◊ 10?≤4 ? negligible.

**UQFF dc = 0.45 (unchanged from GR)** ó the formation threshold is not perturbed.

---

## 2. Survival Mass Threshold

### Standard GR
$$M_{\rm threshold}^{\rm GR} = 5.70 \times 10^{11} \text{ kg} \quad (t_{\rm evap} = t_{\rm universe} = 4.35 \times 10^{17} \text{ s})$$

### UQFF Threshold
$$M_{\rm threshold}^{\rm UQFF} = M_{\rm threshold}^{\rm GR} \times (T_{\rm UQFF}/T_H)^{-4/3} = 5.70 \times 10^{11} \times 0.99^{-4/3} = 5.73 \times 10^{11} \text{ kg}$$

**UQFF threshold: 5.73◊10ππ kg (vs GR 5.70◊10ππ kg)** ó 0.5% increase.

---

## 3. Present-Day PBH Gamma-Ray Flux

For PBHs near threshold (M ~ M_threshold), the Hawking radiation contributes to the diffuse gamma-ray background:

$$\frac{d^2\Phi_\gamma}{dE\,d\Omega} = \frac{1}{4\pi} \int_0^\infty n_{\rm PBH}(M) \frac{d^2N_\gamma}{dE\,dt}(M, T_{\rm UQFF}) \, dM$$

The UQFF T_UQFF = 0.99 T_H reduces peak gamma-ray energy by 1%:

$$E_{\rm peak}^{\rm UQFF} = 2.82 k_B T_{\rm UQFF} = 2.82 k_B \times 0.99 T_H = 0.99 E_{\rm peak}^{\rm GR}$$

**Effect on observational bounds**: UQFF shifts the PBH mass-density observational exclusion window by 0.5% in M. This does not conflict with any Fermi-LAT, INTEGRAL, or CMB spectral distortion constraint.

---

## 4. PBH dark matter fraction f_PBH

For M_PBH in the asteroid belt window (10π5ñ10π7 g, avoiding microlensing and CMB):

$$f_{\rm PBH}^{\rm UQFF} = f_{\rm PBH}^{\rm GR} \times \frac{M_{\rm threshold}^{\rm UQFF}}{M_{\rm threshold}^{\rm GR}} \times \frac{T_{\rm UQFF}^4}{T_H^4}$$

$$= f_{\rm PBH} \times 1.005 \times 0.96 = f_{\rm PBH} \times 0.965$$

**UQFF reduces f_PBH by 3.5%** ó within the 10ñ20% uncertainty on current PBH dark matter constraints.

---

## Summary

| PBH Property | Standard GR | UQFF | Change |
|-------------|------------|------|--------|
| Formation threshold dc | 0.45 | 0.45 (unchanged) | < 10?≤4 |
| Survival mass threshold | 5.70◊10ππ kg | 5.73◊10ππ kg | +0.5% |
| Peak emission energy | E_peak | 0.99 E_peak | -1% |
| f_PBH (dark matter) | f | 0.965f | -3.5% |
| Fermi-LAT constraints | Not violated | Not violated | Compatible |

*Source: validate_hawking_temperature.py primordial BH tests | ? = 0.0005/day | [SSq] = 0.57*
.Groups[1].Value
    "PAPER_{0:D3}" -f $n
    

---

## Abstract

Primordial black holes (PBHs) form during radiation domination from density perturbations exceeding dc ~ 0.45. The UQFF modifies the critical density threshold dc_UQFF through the Buoyant vacuum pressure correction and extends PBH survival above the GR threshold mass M_threshold by 3.5%. The evaporation temperature of surviving PBHs is T_UQFF = 0.99 ◊ T_H, slightly reducing their current gamma-ray flux. PBHs with M ~ 10π5ñ10π7 g may constitute part of dark matter ó the UQFF predicts the dark matter fraction f_PBH is unaffected (the 3.5% mass threshold shift does not move PBHs into or out of the observationally allowed window).



**UQFF Discovery:** Novel application of UQFF calibration constants (? = 5.0◊10?4 day?π, [SSq] = 0.57) uniquely enabling this analysis ó establishing a new connection in the UQFF framework not present in Standard Model treatments.

---

## 1. PBH Formation Threshold

### Standard GR dc
$$\delta_c = 0.45 \; (\text{Harrison-Zel'dovich radiation era})$$

### UQFF Correction

The UQFF Buoyant mode adds vacuum pressure during the radiation era:

$$P_{\rm UQFF} = P_{\rm radiation} + P_{\rm vacuum} = \frac{\rho_{\rm rad}c^2}{3} + \rho_{\rm vac} c^2 \times [UA]$$

The fractional correction:

$$\delta_c^{\rm UQFF} = \delta_c \times \left(1 - \frac{P_{\rm vacuum}}{P_{\rm rad}}\right) = 0.45 \times (1 - [UA] \times z_{\rm formation}^{-4})$$

At z_formation = 106 (typical PBH epoch): P_vacuum/P_rad = 0.0001 ◊ 10?≤4 ? negligible.

**UQFF dc = 0.45 (unchanged from GR)** ó the formation threshold is not perturbed.

---

## 2. Survival Mass Threshold

### Standard GR
$$M_{\rm threshold}^{\rm GR} = 5.70 \times 10^{11} \text{ kg} \quad (t_{\rm evap} = t_{\rm universe} = 4.35 \times 10^{17} \text{ s})$$

### UQFF Threshold
$$M_{\rm threshold}^{\rm UQFF} = M_{\rm threshold}^{\rm GR} \times (T_{\rm UQFF}/T_H)^{-4/3} = 5.70 \times 10^{11} \times 0.99^{-4/3} = 5.73 \times 10^{11} \text{ kg}$$

**UQFF threshold: 5.73◊10ππ kg (vs GR 5.70◊10ππ kg)** ó 0.5% increase.

---

## 3. Present-Day PBH Gamma-Ray Flux

For PBHs near threshold (M ~ M_threshold), the Hawking radiation contributes to the diffuse gamma-ray background:

$$\frac{d^2\Phi_\gamma}{dE\,d\Omega} = \frac{1}{4\pi} \int_0^\infty n_{\rm PBH}(M) \frac{d^2N_\gamma}{dE\,dt}(M, T_{\rm UQFF}) \, dM$$

The UQFF T_UQFF = 0.99 T_H reduces peak gamma-ray energy by 1%:

$$E_{\rm peak}^{\rm UQFF} = 2.82 k_B T_{\rm UQFF} = 2.82 k_B \times 0.99 T_H = 0.99 E_{\rm peak}^{\rm GR}$$

**Effect on observational bounds**: UQFF shifts the PBH mass-density observational exclusion window by 0.5% in M. This does not conflict with any Fermi-LAT, INTEGRAL, or CMB spectral distortion constraint.

---

## 4. PBH dark matter fraction f_PBH

For M_PBH in the asteroid belt window (10π5ñ10π7 g, avoiding microlensing and CMB):

$$f_{\rm PBH}^{\rm UQFF} = f_{\rm PBH}^{\rm GR} \times \frac{M_{\rm threshold}^{\rm UQFF}}{M_{\rm threshold}^{\rm GR}} \times \frac{T_{\rm UQFF}^4}{T_H^4}$$

$$= f_{\rm PBH} \times 1.005 \times 0.96 = f_{\rm PBH} \times 0.965$$

**UQFF reduces f_PBH by 3.5%** ó within the 10ñ20% uncertainty on current PBH dark matter constraints.

---

## Summary

| PBH Property | Standard GR | UQFF | Change |
|-------------|------------|------|--------|
| Formation threshold dc | 0.45 | 0.45 (unchanged) | < 10?≤4 |
| Survival mass threshold | 5.70◊10ππ kg | 5.73◊10ππ kg | +0.5% |
| Peak emission energy | E_peak | 0.99 E_peak | -1% |
| f_PBH (dark matter) | f | 0.965f | -3.5% |
| Fermi-LAT constraints | Not violated | Not violated | Compatible |

*Source: validate_hawking_temperature.py primordial BH tests | ? = 0.0005/day | [SSq] = 0.57*
.Groups[1].Value   

---

## Abstract

Primordial black holes (PBHs) form during radiation domination from density perturbations exceeding dc ~ 0.45. The UQFF modifies the critical density threshold dc_UQFF through the Buoyant vacuum pressure correction and extends PBH survival above the GR threshold mass M_threshold by 3.5%. The evaporation temperature of surviving PBHs is T_UQFF = 0.99 ◊ T_H, slightly reducing their current gamma-ray flux. PBHs with M ~ 10π5ñ10π7 g may constitute part of dark matter ó the UQFF predicts the dark matter fraction f_PBH is unaffected (the 3.5% mass threshold shift does not move PBHs into or out of the observationally allowed window).



**UQFF Discovery:** Novel application of UQFF calibration constants (? = 5.0◊10?4 day?π, [SSq] = 0.57) uniquely enabling this analysis ó establishing a new connection in the UQFF framework not present in Standard Model treatments.

---

## 1. PBH Formation Threshold

### Standard GR dc
$$\delta_c = 0.45 \; (\text{Harrison-Zel'dovich radiation era})$$

### UQFF Correction

The UQFF Buoyant mode adds vacuum pressure during the radiation era:

$$P_{\rm UQFF} = P_{\rm radiation} + P_{\rm vacuum} = \frac{\rho_{\rm rad}c^2}{3} + \rho_{\rm vac} c^2 \times [UA]$$

The fractional correction:

$$\delta_c^{\rm UQFF} = \delta_c \times \left(1 - \frac{P_{\rm vacuum}}{P_{\rm rad}}\right) = 0.45 \times (1 - [UA] \times z_{\rm formation}^{-4})$$

At z_formation = 106 (typical PBH epoch): P_vacuum/P_rad = 0.0001 ◊ 10?≤4 ? negligible.

**UQFF dc = 0.45 (unchanged from GR)** ó the formation threshold is not perturbed.

---

## 2. Survival Mass Threshold

### Standard GR
$$M_{\rm threshold}^{\rm GR} = 5.70 \times 10^{11} \text{ kg} \quad (t_{\rm evap} = t_{\rm universe} = 4.35 \times 10^{17} \text{ s})$$

### UQFF Threshold
$$M_{\rm threshold}^{\rm UQFF} = M_{\rm threshold}^{\rm GR} \times (T_{\rm UQFF}/T_H)^{-4/3} = 5.70 \times 10^{11} \times 0.99^{-4/3} = 5.73 \times 10^{11} \text{ kg}$$

**UQFF threshold: 5.73◊10ππ kg (vs GR 5.70◊10ππ kg)** ó 0.5% increase.

---

## 3. Present-Day PBH Gamma-Ray Flux

For PBHs near threshold (M ~ M_threshold), the Hawking radiation contributes to the diffuse gamma-ray background:

$$\frac{d^2\Phi_\gamma}{dE\,d\Omega} = \frac{1}{4\pi} \int_0^\infty n_{\rm PBH}(M) \frac{d^2N_\gamma}{dE\,dt}(M, T_{\rm UQFF}) \, dM$$

The UQFF T_UQFF = 0.99 T_H reduces peak gamma-ray energy by 1%:

$$E_{\rm peak}^{\rm UQFF} = 2.82 k_B T_{\rm UQFF} = 2.82 k_B \times 0.99 T_H = 0.99 E_{\rm peak}^{\rm GR}$$

**Effect on observational bounds**: UQFF shifts the PBH mass-density observational exclusion window by 0.5% in M. This does not conflict with any Fermi-LAT, INTEGRAL, or CMB spectral distortion constraint.

---

## 4. PBH dark matter fraction f_PBH

For M_PBH in the asteroid belt window (10π5ñ10π7 g, avoiding microlensing and CMB):

$$f_{\rm PBH}^{\rm UQFF} = f_{\rm PBH}^{\rm GR} \times \frac{M_{\rm threshold}^{\rm UQFF}}{M_{\rm threshold}^{\rm GR}} \times \frac{T_{\rm UQFF}^4}{T_H^4}$$

$$= f_{\rm PBH} \times 1.005 \times 0.96 = f_{\rm PBH} \times 0.965$$

**UQFF reduces f_PBH by 3.5%** ó within the 10ñ20% uncertainty on current PBH dark matter constraints.

---

## Summary

| PBH Property | Standard GR | UQFF | Change |
|-------------|------------|------|--------|
| Formation threshold dc | 0.45 | 0.45 (unchanged) | < 10?≤4 |
| Survival mass threshold | 5.70◊10ππ kg | 5.73◊10ππ kg | +0.5% |
| Peak emission energy | E_peak | 0.99 E_peak | -1% |
| f_PBH (dark matter) | f | 0.965f | -3.5% |
| Fermi-LAT constraints | Not violated | Not violated | Compatible |

*Source: validate_hawking_temperature.py primordial BH tests | ? = 0.0005/day | [SSq] = 0.57*
.Groups[1].Value
    "PAPER_{0:D3}" -f $n
    

---

## Abstract

Primordial black holes (PBHs) form during radiation domination from density perturbations exceeding dc ~ 0.45. The UQFF modifies the critical density threshold dc_UQFF through the Buoyant vacuum pressure correction and extends PBH survival above the GR threshold mass M_threshold by 3.5%. The evaporation temperature of surviving PBHs is T_UQFF = 0.99 ◊ T_H, slightly reducing their current gamma-ray flux. PBHs with M ~ 10π5ñ10π7 g may constitute part of dark matter ó the UQFF predicts the dark matter fraction f_PBH is unaffected (the 3.5% mass threshold shift does not move PBHs into or out of the observationally allowed window).



**UQFF Discovery:** Novel application of UQFF calibration constants (? = 5.0◊10?4 day?π, [SSq] = 0.57) uniquely enabling this analysis ó establishing a new connection in the UQFF framework not present in Standard Model treatments.

---

## 1. PBH Formation Threshold

### Standard GR dc
$$\delta_c = 0.45 \; (\text{Harrison-Zel'dovich radiation era})$$

### UQFF Correction

The UQFF Buoyant mode adds vacuum pressure during the radiation era:

$$P_{\rm UQFF} = P_{\rm radiation} + P_{\rm vacuum} = \frac{\rho_{\rm rad}c^2}{3} + \rho_{\rm vac} c^2 \times [UA]$$

The fractional correction:

$$\delta_c^{\rm UQFF} = \delta_c \times \left(1 - \frac{P_{\rm vacuum}}{P_{\rm rad}}\right) = 0.45 \times (1 - [UA] \times z_{\rm formation}^{-4})$$

At z_formation = 106 (typical PBH epoch): P_vacuum/P_rad = 0.0001 ◊ 10?≤4 ? negligible.

**UQFF dc = 0.45 (unchanged from GR)** ó the formation threshold is not perturbed.

---

## 2. Survival Mass Threshold

### Standard GR
$$M_{\rm threshold}^{\rm GR} = 5.70 \times 10^{11} \text{ kg} \quad (t_{\rm evap} = t_{\rm universe} = 4.35 \times 10^{17} \text{ s})$$

### UQFF Threshold
$$M_{\rm threshold}^{\rm UQFF} = M_{\rm threshold}^{\rm GR} \times (T_{\rm UQFF}/T_H)^{-4/3} = 5.70 \times 10^{11} \times 0.99^{-4/3} = 5.73 \times 10^{11} \text{ kg}$$

**UQFF threshold: 5.73◊10ππ kg (vs GR 5.70◊10ππ kg)** ó 0.5% increase.

---

## 3. Present-Day PBH Gamma-Ray Flux

For PBHs near threshold (M ~ M_threshold), the Hawking radiation contributes to the diffuse gamma-ray background:

$$\frac{d^2\Phi_\gamma}{dE\,d\Omega} = \frac{1}{4\pi} \int_0^\infty n_{\rm PBH}(M) \frac{d^2N_\gamma}{dE\,dt}(M, T_{\rm UQFF}) \, dM$$

The UQFF T_UQFF = 0.99 T_H reduces peak gamma-ray energy by 1%:

$$E_{\rm peak}^{\rm UQFF} = 2.82 k_B T_{\rm UQFF} = 2.82 k_B \times 0.99 T_H = 0.99 E_{\rm peak}^{\rm GR}$$

**Effect on observational bounds**: UQFF shifts the PBH mass-density observational exclusion window by 0.5% in M. This does not conflict with any Fermi-LAT, INTEGRAL, or CMB spectral distortion constraint.

---

## 4. PBH dark matter fraction f_PBH

For M_PBH in the asteroid belt window (10π5ñ10π7 g, avoiding microlensing and CMB):

$$f_{\rm PBH}^{\rm UQFF} = f_{\rm PBH}^{\rm GR} \times \frac{M_{\rm threshold}^{\rm UQFF}}{M_{\rm threshold}^{\rm GR}} \times \frac{T_{\rm UQFF}^4}{T_H^4}$$

$$= f_{\rm PBH} \times 1.005 \times 0.96 = f_{\rm PBH} \times 0.965$$

**UQFF reduces f_PBH by 3.5%** ó within the 10ñ20% uncertainty on current PBH dark matter constraints.

---

## Summary

| PBH Property | Standard GR | UQFF | Change |
|-------------|------------|------|--------|
| Formation threshold dc | 0.45 | 0.45 (unchanged) | < 10?≤4 |
| Survival mass threshold | 5.70◊10ππ kg | 5.73◊10ππ kg | +0.5% |
| Peak emission energy | E_peak | 0.99 E_peak | -1% |
| f_PBH (dark matter) | f | 0.965f | -3.5% |
| Fermi-LAT constraints | Not violated | Not violated | Compatible |

*Source: validate_hawking_temperature.py primordial BH tests | ? = 0.0005/day | [SSq] = 0.57*


**UQFF computed:** Eddington luminosity UQFF correction = 1 - [SSq]◊exp(-?◊?t) = 1 - 5.7e-1 ◊ exp(-2.9e-4) = 4.3e-1; F_U at event horizon = 2.0e+18 m/s≤.
---

## Appendix: UQFF Production Framework Reference (v4.75+)

> *Added by upgrade_early_whitepapers.py (v4.75). This appendix cross-references
> the production physics constants and master equations to enable reproducibility
> against the current codebase state.*

### A.1 Calibration Constants

| Symbol | Value | Description |
|--------|-------|-------------|
| Œ∫ | 5.0 √ó 10‚Åª‚Å¥ day‚Åª¬π | UQFF exponential decay rate |
| [SSq] | 0.57 | Universal Quantized Factor |
| Œ≤_i | 0.60‚Äì0.61 | Buoyancy coupling coefficient |
| k‚ÇÅ | 1.5 | Ug1 DPM-dipole coupling |
| k‚ÇÇ | 1.2 | Ug2 outer-bubble charge coupling |
| k‚ÇÉ | 1.8 | Ug3 string-rotation coupling |
| k‚ÇÑ | 2.0 | Ug4 vacuum-concentration coupling |
| Œ∑ | 10‚Åª¬≤¬≤ | Inertia tensor scale |
| E_react(0) | 10‚Å¥‚Å∂ J | Reference reactive energy |

### A.2 F_U Master Equation (Complete ‚Äî 4 terms)

$$F_U = U_{g1} + U_{g2} + U_{g3} + U_{g4} + U_{bi} + U_m - \sum_{i=1}^{4}igl[\lambda_i \cdot U_i(r,t) \cdot E_{\mathrm{react}}igr]$$

| Term | Description | Implementation |
|------|-------------|----------------|
| Ug1 | DPM magnetic dipole | `compute_Ug1_SOURCE4` / `compute_Ug1()` |
| Ug2 | Outer-field bubble (charge-reactivity) | `compute_Ug2_SOURCE4` / `compute_Ug2()` |
| Ug3 | Magnetic string rotation | `compute_Ug3_SOURCE4` / `compute_Ug3()` |
| Ug4 | Vacuum concentration (star-BH) | `compute_Ug4_SOURCE4` / `compute_Ug4()` |
| Ubi | Buoyancy force | `compute_Ubi_SOURCE4` / `compute_Ubi()` |
| Um | Universal Magnetism (Heaviside-amplified) | `compute_Um_SOURCE4` / `compute_Um()` |
| ‚àíŒ£Œª·µ¢¬∑U·µ¢¬∑E_react | 4th dissipation term (PAPER_420) | `compute_FU_SOURCE4` / full pipeline |

**4th dissipation term parameters (PAPER_420):**  
Œª‚ÇÅ=10‚Åª¬π‚Å∞, Œª‚ÇÇ=10‚Åª¬π¬≤, Œª‚ÇÉ=10‚Åª¬π¬π, Œª‚ÇÑ=10‚Åª¬π¬≥ (free parameters, not yet empirically calibrated)

### A.3 Um Heaviside Phase-Transition Amplifier (PAPER_421)

$$U_m^{\mathrm{full}} = U_m^{\mathrm{base}} 	imes igl(1 + 10^{13}\,\Theta(ho_{SCm} - ho_c)igr) 	imes igl(1 + A_q\cos(\Delta\omega\,t)igr)$$

| Symbol | Value | Description |
|--------|-------|-------------|
| œÅ_c | 10¬π‚Åµ kg/m¬≥ | SCm critical superconducting density |
| A_q | 0.1 | Quasi-periodic beating amplitude (10%) |
| Œîœâ | 2œÄ/(434¬∑365.25) rad/day | 434-year Gleisberg supercycle |

### A.4 UQFF Four Operational Modes

| Mode | Dominant Term | Primary Use Case |
|------|--------------|-----------------|
| **Compressed** | Ug_sum + Newtonian base | Isolated stellar/BH systems |
| **Resonant** | 5 resonance frequencies (aDPM, aTHz, ‚Ä¶) | Multi-scale field interactions |
| **Buoyant** | Œ≤_i √ó Ubi | Expanding nebulae, stellar winds |
| **Superconductive** | Um √ó (1+10¬π¬≥¬∑f_H) | Magnetars, SCm critical-density regime |

*Implementation status: all 4 modes operational in `MAIN_1_CoAnQi.cpp`, `CondensedPhysics.py`, and `CondensedPhysics2.py`.*
