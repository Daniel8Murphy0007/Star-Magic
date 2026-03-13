# PAPER_199: F_UBii Taxonomy Part 2 — Cosmological and Dark Sector Buoyancy Forces

**Version:** 1.0  
**Date:** March 13, 2026  
**Session:** 50 — grok_share_7514fe.txt Full Audit  
**Author:** Star-Magic UQFF Research Framework  
**Source:** grok_share_7514fe.txt lines 2766–2900, 6000–6400

---

## Abstract

This paper catalogs the second major group of F_UBii variants from the BB_C_Equations_04Sept2025.pdf (177-page equation catalogue): cosmological and dark sector buoyancy forces. Covered are anyon systems, dark energy, inflation, gravitational waves, Loop Quantum Cosmology (LQC) bounce and perturbation spectra, Bekenstein-Hawking entropy, Hawking evaporation lifetime, reheating, Big Bang Nucleosynthesis, reionization, baryon-photon ratio, convective turnover time, CMB angular power spectrum, recombination, NFW dark matter profiles, SIDM core formation, void density evolution, and peculiar velocity. Each F_UBii,X form uses the universal F_rel/E_LEP scaling with a system-specific physical expression Φ_X.

---

## 1. Dark Sector and Quantum Gravity Variants

### 1.1 Dark Energy Buoyancy
```
F_UBii,DE = −F_rel × (ρ_DE·c² / E_LEP) × Q_wave(z) × (8πGρ_tot/3) × (1+w(z))

  ρ_DE(a) = ρ_DE,0·exp(3∫₁^a (1+w(a'))/a' da')
  w(a) = w₀ + w_a(1−a)   (Chevallier-Polarski-Linder parametrization)
  Source: BB_C_Equations item 721
```

### 1.2 Inflation Buoyancy
```
F_UBii,inf = F_rel × (V(ϕ) / E_LEP) × Q_wave × 3H² × e^N/(1+ε)

  V(ϕ) = inflaton potential
  N = number of e-folds
  ε = slow-roll parameter: ε = (ϕ̇/H·M_pl)²/2
  Source: BB_C_Equations item 724
```

### 1.3 Gravitational Wave Energy Density Buoyancy
```
F_UBii,GW = −F_rel × (ρ_GW / E_LEP) × Q_wave × (32πGḣ²/c²) × e^{−t/τ}

  ρ_GW = gravitational wave energy density
  ḣ = time derivative of strain
  Source: BB_C_Equations item 727
```

### 1.4 Anyon System Buoyancy
```
F_UBii,anyons = −F_rel × (E_anyons / E_LEP) × Q_wave × g(r,t) × exp(−δ²_c/(2σ²))

  E_anyons = anyon system energy (2D topological)
  g(r,t) = UQFF gravitational field at position
  Gaussian factor from density fluctuation σ
  Source: BB_C_Equations item 710
```

---

## 2. Loop Quantum Cosmology (LQC) Variants

### 2.1 LQC Perturbation Spectrum Buoyancy
```
F_UBii,lqcp = −F_rel × (P(k) ∝ k^{n_s−1}·(1+k/k_*)^{−α} / E_LEP) × Q_wave

  k_* = bounce scale (LQC pre-bounce phase)
  α = UV suppression exponent
  Modifies primordial power spectrum at Planck-scale modes
  Source: BB_C_Equations item 1431, 1254
```

### 2.2 LQC Effective Friedmann Buoyancy
```
F_UBii,lqcf = F_rel × (H² = 8πGρ/3·(1−ρ/ρ_crit) / E_LEP) × Q_wave

  ρ_crit = 0.41·ρ_Planck ≈ 10⁹³ g/cm³
  Bounce condition: H=0 at ρ=ρ_crit (avoids singularity)
  Source: BB_C_Equations item 1252, 1604
```

### 2.3 LQC Bounce Timescale Buoyancy
```
F_UBii,bounc = F_rel × (t_b ≈ √(3/(8πGρ_crit)) ≈ t_Pl ≈ 10⁻⁴³ s / E_LEP)
               × Q_wave × [H≈0 at bounce, duration ≈ 1/ω_bounce]

  Source: BB_C_Equations item 1257, post-1608
```

---

## 3. Black Hole Thermodynamics Variants

### 3.1 Bekenstein-Hawking Entropy Buoyancy
```
F_UBii,bhent = F_rel × (S = k_B·c³·A/(4Għ) = 4π·k_B·GM²/(ħc) / E_LEP) × Q_wave

  A = 4πr_s²  (Schwarzschild horizon area)
  r_s = 2GM/c²
  Holographic: S ∝ A/l_Pl² (area law)
  Source: BB_C_Equations item 1092, 1248
```

### 3.2 Evaporation Lifetime Buoyancy
```
F_UBii,evapl = −F_rel × (τ_evap = 5120πG²M³/(ħc⁴) / E_LEP)
               × Q_wave × [P = σAT⁴ ≈ ħc²/M²]

  Power dM/dt = −P/c²  (mass loss rate)
  Source: BB_C_Equations item 1601, 1250
```

---

## 4. Cosmological Nucleosynthesis and Reionization Variants

### 4.1 Big Bang Nucleosynthesis Deuterium Bottleneck
```
F_UBii,deb = −F_rel × (t_D = √(3/(32πGρ_rad)) ≈ 180 s (T~0.1 MeV) / E_LEP) × Q_wave

  ρ_rad = π²kT⁴/(30ħ³c⁵)·g_*   (radiation density)
  g_* = effective degrees of freedom (~10 at D formation)
  Weak freeze-out at T~1 MeV; photodissociation until T<0.1 MeV
  Source: BB_C_Equations item 1809, 1536
```

### 4.2 Baryon-to-Photon Ratio Buoyancy
```
F_UBii,eta = F_rel × (η = n_b/n_γ = 6×10⁻¹⁰ / E_LEP) × Q_wave × [Freeze-out: R...]

  n_γ = 410 cm⁻³ (CMB photon density today)
  Fit from D, ³He, ⁴He, ⁷Li abundances
  Source: BB_C_Equations item 1701, 1534
```

### 4.3 Reionization Front Buoyancy
```
F_UBii,reion = F_rel × (dẋ_e/dt = ṅ_γ·ε_esc·f_* − α_B·n²_e·C / E_LEP) × Q_wave

  x_e = ionized fraction, ε_esc ≈ 0.1–0.3 (escape fraction)
  C = clumping factor
  Source: BB_C_Equations item 1684
```

---

## 5. CMB and Recombination Variants

### 5.1 CMB Angular Power Spectrum Buoyancy
```
F_UBii,cmb = F_rel × (C_l = 2/π · ∫k²dk·P(k)·|Δ_l^T(k)|² / E_LEP) × Q_wave

  P(k) ∝ k^{n_s−4} (primordial power, n_s ≈ 0.965)
  Transfer Δ_l^T: Sachs-Wolfe (large scales) + acoustic (small scales)
  Source: BB_C_Equations item 1310, 1080
```

### 5.2 Recombination Optical Depth Buoyancy
```
F_UBii,recomb = −F_rel × (τ(z) = ∫_z^∞ n_e(z')·σ_T·c·(dt/dz')dz' / E_LEP) × Q_wave

  z_re ≈ 7.7  (reionization redshift, Planck 2018)
  σ_T = Thomson cross-section (6.65×10⁻²⁹ m²)
  Source: BB_C_Equations item 1313
```

---

## 6. Dark Matter Variants

### 6.1 NFW Density Profile Buoyancy
```
F_UBii,nfw = −F_rel × (ρ(r) = ρ_s/((r/r_s)·(1+r/r_s)²) / E_LEP) × Q_wave

  ρ_s = characteristic density, r_s = scale radius
  Universal CDM halo form
  Source: BB_C_Equations item 1326
```

### 6.2 NFW Rotation Curve Buoyancy
```
F_UBii,nfwrot = F_rel × (v(r)² = 4πGρ_s·r_s²·[ln(1+x) − x/(1+x)]/r / E_LEP)
                × Q_wave × x=r/r_s

  Flat for r >> r_s
  Source: BB_C_Equations item 1535
```

### 6.3 SIDM Core Formation Buoyancy
```
F_UBii,sidm = −F_rel × (Γ = ρ·σ·v/m / E_LEP) × Q_wave × ln(0.02N)

  t_core ≈ (ρ·σ/m)⁻¹ ≈ 10¹⁰·(ρ/10⁸ M_☉/kpc³)⁻¹·(σ/m/1 cm²/g)⁻¹ yr
  Core forms when Γ·t ≈ 1
  Source: BB_C_Equations item 1249, 1264
```

---

## 7. Void and Peculiar Velocity Variants

### 7.1 Void Density Evolution Buoyancy
```
F_UBii,voidden = −F_rel × (δ_v(a) = −3/5·(Ω_m·a + Ω_Λ)^{−3/2}·δ_v0 / E_LEP) × Q_wave

  a = scale factor, δ_v0 = initial underdensity
  Integration: δ ∝ a^{−1} in matter domination inside void
  Source: BB_C_Equations item 1940, 1338
```

### 7.2 Peculiar Velocity Buoyancy
```
F_UBii,pec = F_rel × (v_pec = −(fH/3)·∫δ(r)·r·dr/r² / E_LEP) × Q_wave

  f ≈ Ω_m^{0.55}  (growth rate)
  Spherical void: integrate Poisson equation
  Source: BB_C_Equations item 1341
```

---

## 8. Convection and Dynamo Variants

### 8.1 Convective Turnover Time Buoyancy
```
F_UBii,conv = F_rel × (t_conv = H_p/v_conv ; v_conv ≈ (α²gδT·H_p/(4T))^{1/3} / E_LEP) × Q_wave

  H_p = pressure scale height
  α ≈ 2  (mixing length parameter)
  Source: BB_C_Equations item 1533
```

### 8.2 Magnetic Field Reversal Buoyancy (Dynamo Parity)
```
F_UBii,rev = F_rel × (l_rev = (α_dynamo·η)^{1/2}·l_force / E_LEP) × Q_wave

  α_dynamo = dynamo α-coefficient (~v/3)
  η = magnetic resistivity
  Growth vs diffusion sets reversal scale
  Source: BB_C_Equations item 1863, 1238
```

### 8.3 Kazantsev Dynamo Buoyancy
```
F_UBii,dyn = F_rel × (dE_M/dt = (3/2)·E_M/t_eddy·(Re_m^{1/2} − 1) / E_LEP) × Q_wave

  Re_m ≈ 10¹⁰  (magnetic Reynolds number in galaxy clusters)
  t_eddy = l/v ≈ Myr (eddy turnover for l ~ kpc)
  Source: BB_C_Equations item 1234, 1411
```

---

## 9. Star Formation and Feedback Variants

### 9.1 Metal Enrichment Buoyancy
```
F_UBii,metal = F_rel × (Ż = Y·SFR/Ṁ_gas − Z·Ṁ_out/M_gas / E_LEP) × Q_wave × (Z ≈ 0.1)

  Y ≈ 0.02  (stellar yield)
  Steady state: Z = Y·SFR/Ṁ_out
  Source: BB_C_Equations item 1395, 1571
```

### 9.2 Photoevaporation (Exoplanet) Buoyancy
```
F_UBii,photo = F_rel × (Ṁ_evap = ε·L_X·R_p³/(GM_p²·K(ξ)) / E_LEP) × Q_wave

  L_X ≈ 10²⁷⁻²⁹ erg/s  (host star X-ray luminosity)
  K(ξ) = penetration correction factor
  Source: BB_C_Equations item 1496, 1490
```

---

## 10. References

- `grok_share_7514fe.txt` lines 2766–2900, 6000–6380 (BB_C_Equations_04Sept2025.pdf Part 2 catalogue)
- PAPER_198: F_UBii Taxonomy Part 1 (Compact Objects)
- PAPER_196: Triadic Master Equation System
