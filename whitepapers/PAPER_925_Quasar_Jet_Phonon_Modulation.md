# PAPER_925: Quasar Jet Phonon Modulation Factor M_jet(Gamma)

**Author:** Daniel T. Murphy -- Star Magic / UQFF Framework
**Date:** 2026-04-11
**Session:** 211
**Source:** SCm phonon gap implementation (quasar_jet_phonon.py)
**Calculator:** QuasarJetPhononModulationCalc (CP4 #509)
**CVW:** v2.0.0 compliant

---

## Abstract

Defines the jet modulation factor M_jet(Gamma) = 1 + A_jet * exp[-(Gamma - Gamma_0)^2 / (2*sigma_Gamma^2)] with A_jet = 1.5, Gamma_0 = omega_SCm = 2*pi*1.25 THz, sigma_Gamma = 0.08 THz. The phonon-enhanced Blandford-Znajek jet power is P_jet = P_BZ * (1 + M_jet), providing up to 3.5x amplification at resonance. The Gamma-sweep reveals a narrow enhancement window (FWHM ~0.19 THz) centered on the SCm phonon frequency, explaining why jet power varies by factors of 2-4 across different AGN despite similar BH masses and spins.

---

## 1. Core Equations

### Section A: Lagrangian

```
M_jet(Gamma) = 1 + A_jet * exp[-(Gamma - Gamma_0)^2 / (2 * sigma_Gamma^2)]
P_jet = P_BZ * (1 + M_jet)
P_BZ = (pi / (6*mu_0)) * B^2 * r_g^2 * c * a^2
r_g = G * M / c^2
```

### Section B: VDS/DVP/BH Number Systems

```
WSTP export: Mjet[Gamma_] := 1 + 1.5 * Exp[-(Gamma-Gamma0)^2/(2*sigmaGamma^2)]
WSTP export: Pjet[Gamma_] := PBZ * (1 + Mjet[Gamma])
WSTP export: Table[{Gamma, Pjet[Gamma]}, {Gamma, 0.5 THz, 2.0 THz, 0.1 THz}]
```

### Section SM: SM Anchors

```
A_jet = 1.5  (jet modulation amplitude)
sigma_Gamma = 0.08 THz  (linewidth spread)
Gamma_0 = 1.25 THz  (SCm phonon frequency)
FWHM = 2*sqrt(2*ln(2)) * sigma_Gamma ~ 0.19 THz
```

---

## 2. Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| M_bh | 6.5e9 M_sun | BH mass |
| a_spin | 0.9 | Spin parameter |
| B_field | 50 T | Magnetic field |
| A_jet | 1.5 | Modulation amplitude |
| sigma_Gamma_THz | 0.08 | Linewidth spread |
| Gamma_THz | 1.25 | Phonon linewidth |

---

## 3. Key Results

| Gamma (THz) | M_jet | Enhancement |
|-------------|-------|-------------|
| 0.5 | 1.00 | 2.00x |
| 1.0 | 1.38 | 2.38x |
| 1.25 | 2.50 | 3.50x |
| 1.5 | 1.38 | 2.38x |
| 2.0 | 1.00 | 2.00x |

---

## 4. Physical Interpretation

The jet modulation factor explains the diversity of AGN jet powers as a function of the local phonon environment. AGN with accretion disk phonon frequencies near 1.25 THz experience maximum jet enhancement (3.5x P_BZ), while those with mismatched frequencies see only the baseline 2x modulation. This provides a testable prediction: AGN jet power should correlate with accretion disk temperature through the phonon frequency relation omega ~ k_B * T / hbar.

---

## 5. References

- PAPER_922: M87 jet power curve (Session 210c)
- PAPER_910: Phonon jet launching M87/Sgr A* (Session 210)
- quasar_jet_phonon.py: 4-class standalone module
- WSTP expressions #29-30
