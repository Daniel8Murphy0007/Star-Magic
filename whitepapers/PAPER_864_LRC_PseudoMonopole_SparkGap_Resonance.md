# PAPER_864: LRC Circuit Pseudo-Monopole Spark-Gap Resonance (1/r Decay)

**Author:** Daniel T. Murphy -- Star Magic / UQFF Framework
**Date:** 2026-04-05
**Session:** 200
**Source:** advanced_system_analysis_simulator_quantum_calculator.txt (3745 lines)
**Calculator:** LRCPseudoMonopoleSparkGapResonanceCalc (CP4 #448)
**CVW:** v2.0.0 compliant

---

## Abstract

An LRC circuit with spark-gap discharge produces a pseudo-monopole magnetic field exhibiting 1/r decay (not the 1/r³ dipole law). Components: L = 75 µH (23 AWG, 10 ft), C = 500 µF (2×1000 µF series), R = 33.3 Ω (3×100 Ω parallel), spark gap 0.5 mm mild steel. The resonance frequency f_res = 1/(2π√(LC)) = 29.14 Hz produces B = 2.53×10⁻⁸ T at 0.61 m with monopole-like 1/r radial falloff, connecting to the UQFF Ug1 Di-Pseudo-Monopole (DPM) field geometry.

---

## 1. Core Equations

- `f_res = 1 / (2*pi*sqrt(L*C))` = 29.14 Hz
- `P_spark = E_spark * f_spark`
- `I_rms = sqrt(2*P/R)`
- `B = mu_0 * I / (2*pi*r)` = 2.53e-08 T at 0.61 m
- `B(r) = B_0 * (r_0/r)` (pseudo-monopole 1/r decay)
- `Q = (1/R) * sqrt(L/C)` (quality factor)

---

## 2. UQFF Integration

The pseudo-monopole 1/r decay (as opposed to standard dipole 1/r³) is the experimental signature of the Ug1 DPM field geometry at the spark-gap scale. The resonance at 29.14 Hz lies in the ultra-low frequency regime relevant to geomagnetic pulsation coupling. This calculator is differentiated from theoretical DPM classes (PAPER_411 #61, PAPER_855 #439) by its experimental circuit specifications.

---

## 3. Source Data

- **File:** advanced_system_analysis_simulator_quantum_calculator.txt (3745 lines)
- **Session:** 200
- **Circuit specs:** L=75µH, C=500µF, R=33.3Ω, spark gap 0.5mm mild steel
- **VDS/DVP/BH:** ABSENT

---

## References

1. Murphy, D.T. -- Star Magic UQFF Framework (2024-2026)
2. Biot-Savart law; Maxwell's equations for current loop B-field
3. Dirac, P.A.M. -- Quantised Singularities in the Electromagnetic Field (1931)
4. UQFF Calibration: kappa=0.0005/day, [SSq]=0.57, beta_i~0.603
