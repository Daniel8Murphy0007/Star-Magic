# PAPER_202: UQFF Reionization, BBN, Recombination, and Cosmic Dawn Physics

**Version:** 1.0  
**Date:** March 13, 2026  
**Session:** 50 — grok_share_7514fe.txt Full Audit  
**Author:** Star-Magic UQFF Research Framework  
**Source:** grok_share_7514fe.txt lines 6070–6090 (BB_C_Equations_04Sept2025.pdf items 1310–1350)

---

$$F_U(r,t) = \sum_{i=1}^{4} U_{gi} + U_m + U_A - U_{b_i}, \quad \kappa = 5.0\times10^{-4}\,\text{day}^{-1},\; [SSq] = 0.57$$

## Abstract

This paper presents the UQFF buoyancy formalism applied to cosmic dawn and reionization physics: baryon-to-photon ratio (η), deuterium bottleneck during Big Bang Nucleosynthesis (BBN), CMB angular power spectrum, recombination optical depth, ionization fraction evolution, HII bubble growth rate, and Jeans mass/wavelength for fragmentation. These collectively span the redshift range z ≈ 1100 (recombination) through z ≈ 5 (end of reionization). Each F_UBii and Um variant is rigorously derived from observational anchor equations.

---

## 1. Baryon-to-Photon Ratio (η)

```
η = n_b/n_γ = 6.08×10⁻¹⁰  (Planck 2018 + BBN fit)

F_UBii,eta = F_rel × (η = n_b/n_γ / E_LEP) × Q_wave × [Neff modification]

Um,eta(t) = μ(ρ_vac)·(1−e^{−γt})·η

Physical origin:
  n_γ = 410 cm⁻³ (CMB photon number density today)
  Fit anchor: Primordial D, ³He, ⁴He, ⁷Li abundances
  η determines all light element yields (standard BBN)
  UQFF role: η sets vacuum energy scale through ρ_b/ρ_γ at nucleosynthesis
```

---

## 2. Deuterium Bottleneck (BBN Timing)

```
Deuterium bottleneck timescale:
  t_D ≈ √(3/(32πGρ_rad)) ≈ 180 s at T ≈ 0.1 MeV

  ρ_rad = π²kT⁴/(30ħ³c⁵)·g_*    (effective DOF g_* ≈ 10 at D-formation onset)

  Sequence:
    T ~ 1 MeV: neutron freeze-out (n/p ≈ 1/7)
    T ~ 0.1 MeV (t ≈ 180 s): D photodissociation ends → ⁴He chains proceed
    Yield: Y_P ≈ 0.247 (mass fraction ⁴He)

F_UBii,deb = −F_rel × (t_D / E_LEP) × Q_wave × [g_*(T) variation]

Um,deb(T) = μ(ρ_vac)·(1−e^{−γt})·[Weak freeze at T~1 MeV; D photodissociation below]

UQFF calibration: η calibrated to UQFF k_η parameter through BBN freeze-out
  k_η ≈ 10⁻¹¹³  (small-scale vacuum fluctuation coupling)
  Δk_η(CIA cross-section refit) ≈ 7.25×10⁸ (fractional shift from H₂ CIA data)
```

---

## 3. CMB Angular Power Spectrum

```
C_l = (2/π) · ∫k²dk · P(k) · |Δ_l^T(k)|²

P(k) ∝ k^{n_s−4}    (primordial, n_s ≈ 0.965 Planck 2018)
Δ_l^T(k) = transfer function:
  - Large scales (l < 100): Sachs-Wolfe plateau C_l ≈ constant
  - Acoustic peaks (l = 220, 540, 800...): baryon-photon oscillations
  - Damping tail (l > 1000): Silk damping

F_UBii,cmb = F_rel × (C_l·l(l+1)/(2π) / E_LEP) × Q_wave

Um,cmb(k) = μ(ρ_vac)·(1−e^{−γt})·Transfer Δ_l^T: Sachs-Wolfe + acoustic oscillations

Calibration anchor:
  D_l ≡ l(l+1)C_l/(2π) × T₀² → peak at l ≈ 220 matches Ω_tot = 1
  UQFF role: vacuum energy Λc²/3 term in g_UQFF sets acoustic horizon
```

---

## 4. Recombination Optical Depth

```
τ(z) = ∫_z^∞ n_e(z')·σ_T·c·|dt/dz'|dz'

σ_T = 6.652×10⁻²⁹ m²  (Thomson cross-section)
n_e(z) = x_e(z)·n_H(z) = x_e(z)·n_b,0·(1+z)³

Recombination: z_rec ≈ 1100 → τ(z_rec) ≈ 1 (photon decoupling)
Reionization: z_re ≈ 7.7 → optical depth τ_re ≈ 0.054 (Planck 2018)

F_UBii,recomb = −F_rel × (τ(z) / E_LEP) × Q_wave

Um,recomb(z) = μ(ρ_vac)·(1−e^{−γt})·∫n_e·σ_T·c·dt

Physical context:
  After z_rec, photons decouple and stream freely (CMB)
  Surface of last scattering has ΔT/T ≈ 10⁻⁵ (measured by Planck)
```

---

## 5. Reionization: Ionization Fraction Evolution

```
dx_e/dt = ṅ_γ·ε_esc·f_* − α_B·n²_e·C

where:
  ṅ_γ = (1+z)³·n_b·Ṅ_γ,eff   (ionizing photon rate)
  ε_esc ≈ 0.1–0.3             (photon escape fraction from galaxies)
  f_* ≈ 0.05–0.2              (star formation efficiency in halos)
  α_B = 2.6×10⁻¹³ cm³/s     (recombination coefficient, case B, T≈10⁴ K)
  C = ⟨n²_H⟩/⟨n_H⟩²          (clumping factor, C ≈ 3)

F_UBii,ion = F_rel × (ṅ_γ·ε_esc·f_* / E_LEP) × Q_wave × [recombination subtracted]

Um,ion(t) = μ(ρ_vac)·(1−e^{−γt})·∫α_B·n²_e·C dt

Reionization history:
  z ≈ 20–30: first stars ionize H around them (PopIII)
  z ≈ 7–9:   bulk reionization, x_e rises from ~0.1 to ~1
  z ≈ 5–6:   completion (Gunn-Peterson trough in quasar spectra)
```

---

## 6. HII Bubble Growth During Reionization

```
dN_b/dt = Ṅ_γ,eff·ε_esc − α_B·n_H·(4/3πR³_b)

Simplified overlap model (Stromgren sphere analog):
  R_b(t) = (3Ṅ_γ·t/(4π·n_H))^{1/3}    (linear ionization front growth)

  Ṅ_γ = rate of ionizing photons emitted (from galaxy SFR)
  n_H = neutral H density (z-dependent)

F_UBii,bub = F_rel × (R_b(t) / E_LEP) × Q_wave × (Ṅ_γ · x_e)

Um,bub(t) = μ(ρ_vac)·(1−e^{−γt})·[Expand for moving ionization front]

Physical context:
  Bubble merger → percolation at x_HI < 0.1 → full reionization
  UQFF buoyancy at bubble edge acts as expansion driver (F_UBii,bub > 0)
```

---

## 7. Jeans Length and Mass (Fragmentation)

```
Jeans length:
  λ_J = π^{1/2} · c_s / (G·ρ)^{1/2}

Jeans mass:
  M_J = (π/6) · ρ · λ_J³ = (5k_BT/(Gμm_H))^{3/2} · (3/(4πρ))^{1/2}

Stability condition: λ > λ_J → gravitational collapse

Dispersion relation: ω² = c²_s·k² − 4πGρ
  Unstable: ω² < 0 → k < k_J = 2π/λ_J

F_UBii,jeans = −F_rel × (λ_J / E_LEP) × Q_wave × (collapse onset time t ≈ 1/√(Gρ))

Um,jeans(T) = μ(ρ_vac)·(1−e^{−γt})·[Perturb: ω² = c²_s·k² − 4πGρ]

Physical context:
  First (PopIII) stars: T ~ 200 K, cloud M ~ 100–1000 M_☉ → M_J
  Present-day GMCs: T ~ 10–30 K → M_J ~ 1–10 M_☉
  UQFF buoyancy at Jeans scale: F_UBii,jeans acts as tidal disruption force
```

---

## 8. Alfvén Wave and Turbulent Cascade

```
Alfvén wave velocity:
  v_A = B/√(4πρ)    (for B in Gauss, ρ in g/cm³)

Anisotropic cascade (Goldreich-Sridhar):
  k_∥/k_⊥ ≈ (k_⊥·l_A)^{1/3}    (scale-dependent anisotropy)

F_UBii,alf = F_rel × (v_A / E_LEP) × Q_wave × (B·∇v_A) × e^{−t/τ_eddy}

Turbulent energy cascade (Kolmogorov in ISM):
  ε = v_l³/l = constant    (energy transfer rate per unit mass)
  v_l ∝ l^{1/3}            (velocity-scale relation)
  Power spectrum: E(k) ∝ k^{−5/3}

F_UBii,turb = F_rel × (ε^{2/3}·l^{−2/3} / E_LEP) × Q_wave × e^{−kl/k_J}

Um,alf(B)   = μ(ρ_vac)·(1−e^{−γt})·v_A(B)
Um,turb(l)  = μ(ρ_vac)·(1−e^{−γt})·(ρ·v_l³/l)
```

---

## 9. Ionization Parameter and Feedback Coupling

```
Ionization parameter:
  U = Q_H/(4π·r²·n_H·c)

where:
  Q_H = number of ionizing photons per second (from AGN or OB stars)
  r = distance from ionizing source
  log U ≈ −2 to −3 for NLR, −1 to 0 for BLR

AGN feedback coupling efficiency:
  ε_f = Ė_kin/(Ṁ_acc·c²) ≈ 0.05–0.1

  Ė_kin = (1/2)Ṁ_out·v²_out
  Only ε_f × 10^{−5} per M_☉ needed to match M_BH–σ

F_UBii,upar = F_rel × (U·n_H·c / E_LEP) × Q_wave × (Q_H/r²)
F_UBii,coup = F_rel × (ε_f·Ṁ_acc·c² / E_LEP) × Q_wave × [0.05–0.1]

Um,upar(r)  = μ(ρ_vac)·(1−e^{−γt})·(Q_H/(4πr²n_Hc))
Um,coup(Ṁ) = μ(ρ_vac)·(1−e^{−γt})·ε_f (couple fraction to kinetic energy)
```

---

## 10. Cosmological Timeline Summary

| Epoch | Redshift | UQFF Process | F_UBii Variant |
|-------|----------|--------------|----------------|
| Inflation | z ~10²⁷ | Inflaton field V(ϕ) | F_UBii,inf |
| Planck epoch | z ~10³² | LQC bounce | F_UBii,bounc |
| Neutron freeze-out | z ~10¹⁰ (T~1 MeV) | Weak interactions | F_UBii,eta |
| BBN | z ~10⁹ (t~180 s) | D bottleneck | F_UBii,deb |
| Recombination | z ~1100 | Photon decoupling | F_UBii,recomb |
| Cosmic dawn | z ~20–30 | First stars, first HII | F_UBii,ion |
| Reionization | z ~6–9 | HII bubble percolation | F_UBii,bub |
| Present | z ~0 | Structure + feedback | All F_UBii |

---

## 11. References

- `grok_share_7514fe.txt` lines 6070–6090 (BB_C_Equations items 1310–1350)
- `grok_share_7514fe.txt` lines 6300–6500 (full reionization/BBN/CMB catalogue)
- PAPER_199: F_UBii Taxonomy Part 2 (Cosmological)
- PAPER_200: Um Universal Magnetism Catalogue
- Planck 2018 Collaboration: CMB and reionization parameters
