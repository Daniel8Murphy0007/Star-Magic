# PAPER_198: F_UBii Taxonomy Part 1 — Compact Object and Stellar Physics Buoyancy Forces

**Version:** 1.0  
**Date:** March 13, 2026  
**Session:** 50 — grok_share_7514fe.txt Full Audit  
**Author:** Star-Magic UQFF Research Framework  
**Source:** grok_share_7514fe.txt lines 2443–2680 (BB_C_Equations_04Sept2025.pdf catalogue)

---

## Abstract

The UQFF buoyancy force F_UBii is applied across all major astrophysical phenomena by embedding each system's characteristic energy or length scale into the universal F_rel/E_LEP scaling framework. This paper catalogs 18 unique F_UBii variants relating to compact objects and stellar physics: MHD dynamo buoyancy, terminal velocity, Press-Schechter, Hawking radiation, quasi-normal mode ringdown, Blandford-Znajek jet power, Arnett supernova, entanglement, jet velocity, planet migration, AGN feedback, angular momentum accretion, J-type shock, Sedov-Taylor blast wave, GRB afterglow, SIDM core formation, ionization fronts, superfluid glitch, and virial theorem.

---

## 1. General F_UBii Framework

```
F_UBii,X = ±F_rel × (Φ_X / E_LEP) × Q_wave × [decay/oscillation factor]

where:
  F_rel ≈ 4.3×10³³ N  (relativistic UQFF force scale)
  E_LEP = energy of characteristic LEP scale (J)
  Q_wave = quantum wave amplitude (J/m³), std ≈ 6.33×10⁴
  Φ_X = system-specific physical expression
```

---

## 2. Compact Object Variants

### 2.1 MHD Dynamo Buoyancy
```
F_UBii,MHD = F_rel × (E_M/t_eddy/E_LEP) × Q_wave × (3/2)(Re_m^{1/2} − 1)

  E_M = magnetic energy density
  t_eddy = l/v ~ Myr (eddy turnover time)
  Re_m ≈ 10¹⁰    (magnetic Reynolds number)
  Source: BB_C_Equations item 748
```

### 2.2 Terminal Velocity Buoyancy
```
F_UBii,termv = F_rel × (√(2GM(1−Γ)/r_launch) / E_LEP) × Q_wave × (τ·L/c)

  Γ = L/L_Edd ≈ 1   (Eddington ratio for wind-driven systems)
  r_launch ≈ 100 R_s  (wind launch radius)
  τ = optical depth
  Source: BB_C_Equations item 853, 1823
```

### 2.3 Hawking Radiation Buoyancy
```
F_UBii,haw = F_rel × (ħc³/(8πGMk_B) / E_LEP) × Q_wave × (κ/(2π))

  T_H = ħc³/(8πGMk_B)  (Hawking temperature)
  κ = c⁴/(4GM)          (surface gravity)
  Source: BB_C_Equations item 901, 1246
```

### 2.4 Quasi-Normal Mode Ringdown Buoyancy
```
F_UBii,qnm = −F_rel × (c³/(2πGM) · (0.3737 + 0.088·a_f) / E_LEP)
              × Q_wave × e^{−t/τ} × [τ ∝ M]

  a_f ≈ 0.69   (dimensionless BH spin post-merger)
  τ = ringdown decay time
  f_QNM = c³/(2πGM_f) · f(0.3737 + 0.088·a_f) (Berti fits l=2, m=2 mode)
  Source: BB_C_Equations item 945, 1293
```

### 2.5 Blandford-Znajek Jet Power Buoyancy
```
F_UBii,bz = F_rel × (1/32·B²R_H⁴Ω_H²/c / E_LEP)
             × Q_wave × (ac/(2R_H)) × (1+t_var)

  R_H = GM/c²  (horizon radius)
  Ω_H = ac/(2R_H)  (BH angular velocity)
  Source: BB_C_Equations item 1147

F_UBii,bz2 = F_rel × (κ/(16π) · Φ²_BH · Ω²_BH/c / E_LEP)
              × Q_wave × (ac³/(2GM)) × 0.05−1

  Φ_BH = B·π·r_H²  (BH magnetic flux thread; EHT-calibrated)
  Source: BB_C_Equations item 967, 1316
```

### 2.6 Arnett Supernova Buoyancy
```
F_UBii,arnett = F_rel × (M_Ni·ε_Ni/t_d / E_LEP) × Q_wave × e^{−t/τ}

  M_Ni = nickel mass (~0.5 M_☉)
  ε_Ni = nickel decay energy
  t_d² = 3κM/(4πcv²)  (diffusion time)
  Source: BB_C_Equations item 1035
```

### 2.7 TOV Hydrostatic Equilibrium Buoyancy
```
F_UBii,tov = −F_rel × (dP/dr = −Gm(r)ρ(r)/r² · (1 + P(r)/(ρ(r)c²))
              × (1 + 4πr³P(r)/(m(r)c²)) / (1 − 2Gm(r)/(rc²)) / E_LEP) × Q_wave

  Includes all GR corrections (Schwarzschild metric, pressure-energy coupling)
  Source: BB_C_Equations item 1300
```

### 2.8 Pulsar Spin-Down Buoyancy
```
F_UBii,puls = −F_rel × (τ = P/(2Ṗ) / E_LEP) × Q_wave

  P = period, Ṗ ≈ 10⁻¹⁵ s/s (period derivative, Vela)
  I·dΩ/dt = −L/Ω   (torque equation)
  Source: BB_C_Equations item 1302
```

---

## 3. Stellar/Accretion Variants

### 3.1 Jet Velocity Buoyancy
```
F_UBii,jetvel = F_rel × (v_j ≈ v_K·(r_A/r_0)^{1/2} / E_LEP)
                × Q_wave × (B/√(4πρ)) × (1+t/τ_A)

  v_K = √(GM/r_0)  (Keplerian velocity at footpoint, r_0 = 1–10 AU)
  r_A = Alfvén radius (10–50 AU, from POETS protostellar data)
  Source: BB_C_Equations item 1096, 1272
```

### 3.2 Type-I Planet Migration Buoyancy
```
F_UBii,migr = −F_rel × (−f·2(GM_p)² · Σ/(M_*·Ω·a⁵·(H/r)³) / E_LEP)
               × Q_wave × [τ = M_p/Ṁ_acc ≈ 10⁶ yr]

  f ≈ −1.36  (Lindblad/corotation factor, inward migration)
  Σ = disk surface density
  Source: BB_C_Equations item 1121, 1322
```

### 3.3 Superfluid Vortex Glitch Buoyancy
```
F_UBii,glitch = F_rel × (Δν = I_s/I · ν₀ · (1−e^{−t/τ_q}) / E_LEP)
                × Q_wave × ΔΩ × e^{−t/τ_q}

  I_s = superfluid moment of inertia (~10³⁸ kg·m²)
  τ_q = quench timescale
  ΔΩ = angular velocity jump
  Source: BB_C_Equations item 1753, 1304
```

---

## 4. Shock and Blast Wave Variants

### 4.1 J-Type Shock Buoyancy (Rankine-Hugoniot)
```
F_UBii,jshock = F_rel × ((ρ₁v₁² + P₁) / E_LEP)
                × Q_wave × (v₂/v₁) × (γ+1)/(γ−1+2/M²)

  γ = 5/3  (adiabatic index)
  M = v₁/c_s  (Mach number, from Chandra X-rays in HH 154)
  Source: BB_C_Equations item 1193, 1274
```

### 4.2 Sedov-Taylor Blast Wave Buoyancy
```
F_UBii,sedov = F_rel × ((E·t²/ρ)^{1/5} / E_LEP) × Q_wave × [d/dt(Mv)=0] × t^{2/5}

  E ≈ 10⁵¹ erg  (explosion energy)
  ρ = ambient density
  R(t) = (Et²/ρ)^{1/5}  (self-similar blast radius)
  Source: BB_C_Equations item 1207, 1288
```

### 4.3 C-Type Shock Buoyancy (Magnetized)
```
F_UBii,cshock = −F_rel × ((∇×B)×B/(4π) + ρ_i·ν_in·(v_n−v_i) / E_LEP) × Q_wave

  C-shocks: continuous, multi-fluid MHD (ions+neutrals+magnetic field coupled)
  ν_in = ion-neutral collision frequency
  Source: BB_C_Equations item 1276
```

---

## 5. Feedback and Outflow Variants

### 5.1 AGN Feedback Buoyancy
```
F_UBii,agn = F_rel × (f(v_out)·L_AGN/c / E_LEP) × Q_wave × (Ṁ_out·v_out) × v_out⁻¹

  Momentum-driven: p_term = Ṁ_out·v_out = f(v_out)·L_AGN/c
  Ṁ_out ≈ 10–100 M_☉/yr
  Source: BB_C_Equations item 1165, 1314
```

### 5.2 Angular Momentum Accretion Buoyancy
```
F_UBii,ang = −F_rel × (Ṁ·r²·Ω / E_LEP) × Q_wave × T_B × e^{−t/τ_disk}

  T_B = B²r³/(4π)  (magnetic braking torque)
  τ_disk = disk decay timescale
  Source: BB_C_Equations item 1189, 1270
```

### 5.3 Feedback Energy Coupling Buoyancy
```
F_UBii,coup = F_rel × ((1/2)·Ṁ_w·v_w² / (Ṁ_acc·c²·10) / E_LEP) × Q_wave × 0.05–0.1

  ε_f = Ė_kin/(Ṁ_acc·c²) ≈ 0.05–0.1  (coupling fraction)
  Source: BB_C_Equations item 1331, 1554
```

---

## 6. Structure Formation Variants

### 6.1 Press-Schechter Halo Mass Function Buoyancy
```
F_UBii,ps = F_rel × (dn/dM = √(2/π)·(ρ₀/M)·(δ_c/σ(M))·|d ln σ/d ln M|·exp(−δ²_c/(2σ²)) / E_LEP)
             × Q_wave × ΔΩ (Gaussian part from σ ≈ δ_c)

  δ_c ≈ 1.69  (critical collapse overdensity)
  Source: BB_C_Equations item 877, 1574
```

### 6.2 Structure Growth Rate Buoyancy
```
F_UBii,grow = −F_rel × (δ̈ + 2H·δ̇ = (3/2)·Ω_m·H²·δ/a³ / E_LEP) × Q_wave

  D(a) = 5Ω_m/2 · ∫₀^a da'/(a'H(a')/H₀)³  (growth factor)
  Growing mode: D ∝ a in matter era, suppressed by dark energy
  Source: BB_C_Equations item 1371, 1244
```

---

## 7. GRB Variants

### 7.1 GRB Fireball Expansion Buoyancy
```
F_UBii,fire = F_rel × (Γ(r) = r/R₀ (r<R_s), Γ=η (r>R_s) / E_LEP) × Q_wave

  R₀ ≈ 10⁷ cm  (initial fireball radius)
  R_s = η²R₀   (saturation radius)
  Source: BB_C_Equations item 1482, 1306
```

### 7.2 GRB Afterglow Synchrotron Buoyancy
```
F_UBii,after = −F_rel × (F_ν ∝ ν^{−(p−1)/2}·t^{−3(p−1)/4} (ν_m<ν<ν_c) / E_LEP) × Q_wave

  p ≈ 2.2–2.5  (electron power-law index)
  Electrons accelerated by DSA, spectrum N(E) ∝ E^{−p}
  Source: BB_C_Equations item 1227, 1308
```

---

## 8. References

- `grok_share_7514fe.txt` lines 2443–2680 (BB_C_Equations_04Sept2025.pdf Part 1 catalogue)
- PAPER_196: Triadic Master Equation System
- PAPER_197: F_U_Bi_i Extended Integral
