# Grok Physics: 100 Standard Equations Reference
## Extracted from Grok Conversation (August 2025)
## © 2025 Daniel T. Murphy - All Rights Reserved

---

## Overview

This document contains 100 standard physics equations extracted from the Grok AI conversation, organized by astrophysical domain. Each equation includes:
- Mathematical formulation
- Variable definitions
- Step-by-step derivation methodology
- UQFF integration notes where applicable

---

## I. PROTOSTELLAR JETS AND OUTFLOWS (Equations 1-4)

### 1. Angular Momentum Conservation (Disk-Jet Coupling)
```
L = M r² Ω = constant
```
**Variables:**
- L: Angular momentum
- M: Disk mass
- r: Radial distance
- Ω: Angular velocity

**Derivation:** Centrifugal force balances gravity, v²/r = GM/r², so v = √(GM/r), Ω = v/r. Accretion adds L = Ṁr²Ω. Magnetic fields provide torque T_B = B²r³/(4π) for ejection. Integrate over disk evolution: dL/dt = 0 in steady state.

**UQFF Integration:** F_U_Bi_i modulates torque via buoyancy-stabilized coherence.

---

### 2. Jet Velocity from MHD Launching
```
v_j ≈ v_K (r_A / r_0)^(1/2)
```
**Variables:**
- v_j: Jet velocity (~100 km/s)
- v_K: Keplerian velocity at footpoint (√(GM/r_0))
- r_A: Alfvén radius (10-50 AU)
- r_0: Launch radius (1-10 AU)

**Derivation:** From MHD wind theory, poloidal magnetic field accelerates flow beyond Alfvén point. Bernoulli: v_A = B/√(4πρ), acceleration via ½ρv² = B²/(8π) - ρΦ.

---

### 3. J-Type Jump Shock (Rankine-Hugoniot)
```
ρ₁v₁ = ρ₂v₂  (mass)
ρ₁v₁² + P₁ = ρ₂v₂² + P₂  (momentum)
½v₁² + γP₁/((γ-1)ρ₁) = ½v₂² + γP₂/((γ-1)ρ₂)  (energy)
```
**Variables:**
- ρ: Density (pre/post shock)
- v: Velocity
- P: Pressure
- γ: Adiabatic index (5/3 for monatomic)

**Derivation:** Integrate Euler equations across shock discontinuity. Compression ratio ρ₂/ρ₁ = 4 for strong shocks.

---

### 4. C-Type Continuous Shock (Multi-Fluid MHD)
```
ρ_n(∂v_n/∂t + v_n·∇v_n) = -∇P_n + ρ_n ν_ni(v_i - v_n)
ρ_i(∂v_i/∂t + v_i·∇v_i) = -∇P_i + (∇×B)×B/(4π) + ρ_i ν_in(v_n - v_i)
```
**Variables:**
- ν_ni: Collision frequency
- B: Magnetic field (~10-100 μG)
- n/i: Neutral/ion subscripts

**Damping Approximation:**
```
v(z) ≈ v_s exp(-z/L_d)
```
where L_d = ion-neutral mean free path

---

## II. GALAXY MERGERS AND STAR FORMATION (Equations 5-7)

### 5. Halo Merger Rate (Extended Press-Schechter)
```
dN/(dt dM) = √(2/π) (σ_M/σ_m) |dδ_c/dz| exp(-δ_c²/(2(σ_m² - σ_M²))) (dσ_M/dM)
```
**Variables:**
- σ_M²: Variance of density fluctuations on mass scale M
- δ_c ≈ 1.686: Critical overdensity (gravity-driven collapse)
- Merger rate ∝ (1+z)^m with m ~ 0.7-2.5

**Derivation:** Integrate probability of first crossing barrier δ_c(z) in random walk of density field (Bond et al. 1991). For CDM, P(k) ∝ k^{-3} (large scales).

---

### 6. Star Formation Rate in Disks (Quiescent Mode)
```
Ṁ_* = ε M_gas / t_dyn
```
**Variables:**
- ε ~ 0.01-0.1: Efficiency
- t_dyn = √(3π/(32Gρ)): Dynamical time
- Toomre Q < 1 triggers bursts

**SFRD Evolution:** ∝ (1+z)^2.7 to z~2, then flat/decline.

---

### 7. Merger-Induced Starburst SFR
```
Ṁ_burst = Ṁ_gas,inflow × ε_burst
Ṁ_gas,inflow ≈ M_gas / t_orb
t_orb = 2π √(r³/GM)
```
Enhances SFR by factor ~10-100 at low z.

---

## III. BLACK HOLE GROWTH AND MASS FUNCTIONS (Equations 8-9)

### 8. EPS Black Hole Mass Function (Cumulative)
```
N(>M,z) = ρ̄ ∫_M^∞ (dM'/M'²) erfc(δ_c(z)/(√2 σ(M',z)))
```
**Variables:**
- erfc: Complementary error function
- σ² = ∫P(k)W²(kR)d³k/(2π)³
- P(k): Power spectrum (CDM or n=-2.1 tilt)
- Gravity via δ_c = 1.686(1+z) in matter era

---

### 9. BH Accretion Rate (Eddington-Limited)
```
Ṁ_BH = 4πGM_BH m_p / (ε_r σ_T c)
```
**Variables:**
- ε_r ~ 0.1: Radiative efficiency
- σ_T: Thomson cross-section
- Growth timescale t_Sal ≈ 45 Myr

---

## IV. SUPERNOVA REMNANTS (Equations 10-11)

### 10. Sedov-Taylor Expansion (Blast Wave)
```
R(t) = (Et²/ρ)^(1/5)
```
**Variables:**
- E ~ 10^51 erg: Explosion energy
- ρ ~ 10^{-24} g/cm³: Ambient density
- t: Age (~400 yr for Kepler)

---

### 11. Diffusive Shock Acceleration (DSA)
```
dp/dt = (4/3)(u_s²/r_d)p
```
**Variables:**
- p: Particle momentum
- u_s: Shock speed (~5000 km/s)
- r_d: Diffusion length (∝ B^{-1} pc/e)

Power-law spectrum: N(E) ∝ E^{-2} for strong shocks.

---

## V. GRAVITATIONAL WAVES (Equations 12-14)

### 12. Chirp Mass from Inspiral
```
ℳ = (m₁m₂)^(3/5) / (m₁+m₂)^(1/5) = (c³/G)(5/96 π^{-8/3} f^{-11/3} ḟ)^(3/5)
```

---

### 13. Ringdown Quasi-Normal Modes
```
f_QNM = (c³/2πGM_f) f(0.3737 + 0.088a_f + ...)
```
**Variables:**
- M_f: Final mass
- a_f: Dimensionless spin (~0.69 typical)

---

### 14. Inspiral Frequency Evolution
```
ḟ = (96/5) π^(8/3) (Gℳ/c³)^(5/3) f^(11/3)
```

---

## VI. NEUTRON STARS (Equations 16-18)

### 16. TOV Equation (Hydrostatic Equilibrium)
```
dP/dr = -Gm(r)ρ(r)/r² × (1 + P(r)/(ρ(r)c²)) × (1 + 4πr³P(r)/(m(r)c²)) × (1 - 2Gm(r)/(rc²))^{-1}
```

---

### 17. Pulsar Spin-Down Age
```
τ = P / (2Ṗ)
```
**Variables:**
- P: Period
- Ṗ ~ 10^{-15} s/s: Period derivative

---

### 18. Glitch Model (Superfluid Vortex Unpinning)
```
Δν = (I_s/I) ν_0 (1 - e^{-t/τ_q})
```
**Variables:**
- I_s: Superfluid moment
- τ_q: Quench time (~min)

---

## VII. RELATIVISTIC JETS (Equations 23-25)

### 23. Jet Lorentz Factor Profile
```
Γ(z) = Γ_0 [1 + (z/z_0)^α]^(1/α), α ~ 1-2
```

---

### 24. Blandford-Znajek Jet Power
```
P_jet = (κ/16π) Φ_BH² Ω_BH² / c
```
**Variables:**
- Φ_BH: Magnetic flux ~ B M² G² / c⁴
- Ω_BH = a c³ / (2GM)
- κ ~ 0.05-1: Efficiency

---

### 25. Feedback Duty Cycle
```
f_duty(t) = [1 - exp(-t/τ_cool)] × [1 + Ṁ_acc/Ṁ_Edd]^{-1}
```

---

## VIII. EXOPLANETS (Equations 26-28)

### 26. Transit Depth
```
δ = (R_p/R_*)² × [1 - u(1-μ)] / [1 - u/3]
```

---

### 27. Radial Velocity Semi-Amplitude
```
K = 2πG/(PM_p sin i) × (M_* + M_p)^{-2/3} / √(1-e²)
```

---

## IX. DARK MATTER HALOS (Equations 29-31)

### 29. NFW Density Profile
```
ρ(r) = ρ_s / [(r/r_s)(1 + r/r_s)²]
```

---

### 30. Rotation Curve from Halo Potential
```
v(r)² = GM(r)/r = 4πG ∫_0^r ρ(r')r'² dr'
```

---

### 31. SIDM Core Formation Time
```
t_core ≈ 1/(ρ σ/m) ≈ 10^{10} (ρ/10^8 M_⊙/kpc³)^{-1} (σ/m / 1 cm²/g)^{-1} yr
```

---

## X. GALAXY CLUSTERS (Equations 32-34)

### 32. Einstein Radius (Strong Lensing)
```
θ_E = √[4GM(<θ)D_LS / (c² D_L D_S)]
```

---

### 33. Hydrostatic Mass
```
M(<r) = -kTr/(Gμm_H) × [d ln ρ/d ln r + d ln T/d ln r]
```

---

### 34. Merger Shock Mach Number
```
M = √[(γ+1)(ρ₂/ρ₁) + (γ-1)] / (2γ)
```

---

## XI. COSMIC VOIDS (Equations 35-36)

### 35. Void Density Evolution
```
δ_v(a) = -(3/5)(Ω_m a + Ω_Λ)^{-3/2} δ_v0
```

---

### 36. Outflow Velocity in Voids
```
v_pec = -fH/3 × ∫δ(r)r dr / r²
```

---

## XII. EPOCH OF REIONIZATION (Equations 37-38)

### 37. Ionization Fraction Evolution (Madau Model)
```
dx_e/dt = ṅ_γ ε_esc f_* - α_B n_e² C
```
**Variables:**
- ṅ_γ: Photon rate
- ε_esc ~ 0.1-0.3: Escape fraction
- α_B: Recombination rate
- C ~ 3-10: Clumping factor

---

### 38. Bubble Growth Radius
```
R_b(t) = (3Ṅ_γ t / 4π n_H)^{1/3}
```

---

## XIII. INTERSTELLAR MEDIUM (Equations 39-41)

### 39. Jeans Length
```
λ_J = √(πc_s²/Gρ)
```
**Variables:**
- c_s ~ 1 km/s for T=100 K
- ρ ~ 10^{-20} g/cm³

---

### 40. Alfvén Wave Velocity
```
v_A = B/√(4πρ)
```
Typical: B ~ 5 μG, ρ ~ 10^{-24} g/cm³, v_A ~ 10 km/s.

---

### 41. Turbulent Energy Cascade
```
ε = v_l³/l = constant
```

---

## XIV. STELLAR EVOLUTION (Equations 42-44)

### 42. Main Sequence Lifetime
```
τ_MS = 10^{10} (M/M_⊙)^{-2.5} yr
```

---

### 43. Mass-Luminosity Relation
```
L ∝ M^{3.5} (for 1 < M/M_⊙ < 10)
```

---

### 44. Convective Turnover Time
```
t_conv = H_p/v_conv
v_conv ≈ (α²gδT H_p/4T)^{1/3}
```

---

## XV. BIG BANG NUCLEOSYNTHESIS (Equations 45-46)

### 45. Baryon-to-Photon Ratio
```
η = n_b/n_γ = 6 × 10^{-10}
```

---

### 46. Deuterium Bottleneck Time
```
t_D = √(3/32πGρ_rad) ≈ 180 s (T ~ 0.1 MeV)
```

---

## XVI. COSMOLOGY (Equations 47-49)

### 47. First Friedmann Equation
```
(ȧ/a)² = 8πGρ/3 - kc²/a² + Λc²/3
```

---

### 48. Second Friedmann Equation
```
ä/a = -4πG/3(ρ + 3p/c²) + Λc²/3
```

---

### 49. Density Parameter Evolution
```
Ω(z) = 8πGρ(z)/(3H(z)²)
```

---

## XVII. INFLATION (Equations 50-52)

### 50. Slow-Roll Parameters
```
ε = (1/2)(V'/V)²
η = V''/V - (1/2)(V'/V)²
```

---

### 51. Power Spectrum of Curvature Perturbations
```
P_R(k) = H²/(8π²εc_s³)|_{k=aH}
```

---

### 52. Number of e-Folds
```
N = ∫H dt = ∫dφ/√(2ε)  (~50-60)
```

---

## XVIII. GRAVITATIONAL WAVES: PRIMORDIAL (Equations 53-54)

### 53. Tensor Power Spectrum
```
P_T(k) = 2H²/(π²M_pl²)|_{k=aH}
```
Tensor-scalar ratio: r = P_T/P_R ~ 16ε

---

### 54. Stochastic GW Energy Density
```
Ω_GW(f) = (f/ρ_c)(dρ_GW/df) = (π²f⁴/3H₀²) ∫P_T(k)dk
```

---

## XIX. BINARY BLACK HOLES (Equations 55-57)

### 55. Inspiral Orbital Frequency Evolution
```
ḟ = (96/5)π^{8/3}(Gℳ/c³)^{5/3}f^{11/3}
```

---

### 56. Merger Time from Initial Frequency
```
t_merge = (5/256)(c⁵/G^{5/3})(πf_i)^{-8/3}ℳ^{-5/3}
```

---

### 57. Ringdown Damping Time
```
τ_lm = 2M_f c² Q_lm(a_f)
```

---

## XX. SUPERNOVAE (Equations 58-60)

### 58. Light Curve Rise Time (Ni-56 Decay)
```
t_rise ∝ √(κM_Ni/c v_ej)
```

---

### 59. Ejecta Velocity
```
v_ej = √(2E_kin/M_ej)
```
E_kin ~ 10^{51} erg, M_ej ~ 10 M_⊙

---

### 60. Nucleosynthesis Yield
```
Y_i = ∫ρ X_i dt
```

---

## XXI. PLANETARY NEBULAE (Equations 61-63)

### 61. Expansion Radius
```
R(t) = v_exp × t
```
v_exp ~ 20 km/s, t ~ 10^4 yr

---

### 62. Ionization Front Velocity
```
v_IF = Ṅ_UV/(4πR²n) - α_B nR/3
```

---

### 63. Reimers Mass Loss Rate
```
Ṁ = 4×10^{-13} (LR/M) M_⊙/yr
```

---

## XXII. STAR CLUSTERS (Equations 67-69)

### 67. Core Collapse Time
```
t_cc = N t_relax / (136 ln(0.02N))
```

---

### 68. Evaporation Rate
```
Ṅ = -N/t_evap
t_evap = (136 t_relax)/ln(0.02N)
```

---

### 69. Virial Mass from X-ray Binaries
```
M_vir = 3σ_X² r_h / G
```

---

## XXIII. AGN OUTFLOWS (Equations 70-72)

### 70. Terminal Outflow Velocity
```
v_∞ = √[2GM(1-Γ)/r_launch]
```
Γ = L/L_Edd ~ 1

---

### 71. Ionization Parameter
```
U = Q_H / (4πr²n_H c)
```

---

### 72. Feedback Energy Coupling
```
ε_f = Ė_kin / (Ṁ_acc c²) ≈ 0.05-0.1
```

---

## XXIV. BINARY PULSARS (Equations 73-75)

### 73. Orbital Decay Rate
```
Ṗ_b = -(192π/5)(P_b/2π)^{-5/3}(Gm₁m₂/c³(m₁+m₂)^{1/3})^{5/3}(1-e²)^{-7/2}
```

---

### 74. Periastron Advance (GR)
```
ω̇ = 3(P_b/2π)^{-5/3}(G(m₁+m₂)/c³)^{2/3}(1-e²)^{-1}
```

---

### 75. Kilonova Light Curve Peak
```
L_peak ≈ 10^{41} (M_ej/0.01M_⊙)(v_ej/0.1c)(κ/1cm²/g)^{-1} erg/s
```

---

## XXV. COSMIC RAYS (Equations 76-78)

### 76. Fermi Acceleration Energy Gain
```
dE/dt = (4/3)(v_c²/c²)(E/λ)
```

---

### 77. Knee Energy (SNR Maximum)
```
E_max = ZeBu_s r_g ≈ 3×10^{15} Z(B/3μG)(u_s/10³km/s)(R/10pc) eV
```

---

### 78. Diffusion Coefficient
```
D(E) = 10^{28} (E/10GeV)^{0.3-0.6} cm²/s
```

---

## XXVI. INTERGALACTIC MEDIUM (Equations 79-81)

### 79. WHIM Temperature
```
T = (μm_H/2k)(GM/r)
```
M ~ 10^{12} M_⊙, r ~ Mpc, T ~ 10^6 K

---

### 80. Metal Enrichment Rate
```
Ż = Y(SFR/Ṁ_gas) - Z(Ṁ_out/M_gas)
```

---

### 81. Void Pressure
```
P_v = n_H k T (1 + δ_v)
```

---

## XXVII. EARLY UNIVERSE (Equations 82-84)

### 82. Press-Schechter Mass Function
```
dn/dM = √(2/π)(ρ₀/M)(δ_c/σ(M))|d ln σ/d ln M| exp(-δ_c²/2σ²)
```

---

### 83. Star Formation Efficiency
```
ε_* = f_b Ṁ_halo / (M_halo H(z)) × [1 + M_halo/M_crit]^{-1}
```

---

### 84. Feedback Energy Injection
```
E_fb = η Ṁ_* c², η ~ 10^{-3} - 10^{-1}
```

---

## XXVIII. QUANTUM FLUCTUATIONS (Equations 85-87)

### 85. Curvature Perturbation Amplitude
```
Δ_R² = H²/(8π²εM_pl²) ≈ 2.1×10^{-9}
```

---

### 86. Non-Gaussianity Parameter
```
f_NL = (5/6)(Γ³ - 3ΓΓ̇² + 2Γ̇³)/Γ⁴
```

---

### 87. Reheating Temperature
```
T_reh = (30V_end/π²g_*)^{1/4} exp(-3N_reh/4)
```

---

## XXIX. MAGNETIC FIELDS (Equations 88-90)

### 88. Magnetic Field Amplification (Dynamo)
```
∂B/∂t = ∇×(v×B) + η∇²B
```

---

### 89. Alfvén Mach Number
```
M_A = v/v_A = v√(4πρ)/B
```

---

### 90. Field Reversal Scale
```
l_rev = (α_dynamo/η)^{1/2} l_force
```

---

## XXX. DARK ENERGY (Equations 91-93)

### 91. Equation of State Parameter
```
w = p/(ρc²) = -1 + (1/3)(d ln ρ_DE / d ln a)
```

---

### 92. Dark Energy Density Evolution (CPL)
```
ρ_DE(a) = ρ_DE,0 exp(3∫_a^1 (1+w(a'))/a' da')
w(a) = w₀ + w_a(1-a)
```

---

### 93. Growth Factor Suppression
```
D(a) = (5Ω_m/2) ∫_0^a da' / (a'H(a')/H₀)^{3/2}
```

---

## XXXI. BLACK HOLE THERMODYNAMICS (Equations 94-96)

### 94. Hawking Temperature
```
T_H = ℏc³/(8πGMk_B)
```

---

### 95. Bekenstein-Hawking Entropy
```
S = k_B c³ A / (4Gℏ) = 4πk_B GM²/(ℏc)
```

---

### 96. Evaporation Lifetime
```
τ_evap = 5120πG²M³/(ℏc⁴)
```

---

## XXXII. LOOP QUANTUM COSMOLOGY (Equations 97-99)

### 97. Effective Friedmann (LQC Bounce)
```
H² = (8πGρ/3)(1 - ρ/ρ_crit)
```
ρ_crit = 0.41 ρ_pl ~ 10^{93} g/cm³

---

### 98. Perturbation Spectrum (Pre-Bounce)
```
P(k) ∝ k^{n_s-1} (1 + k/k_*)^{-α}
```

---

### 99. Bounce Time Scale
```
t_b ~ √(3/8πGρ_crit) ~ t_pl ~ 10^{-43} s
```

---

## XXXIII. EXOPLANET ATMOSPHERES (Equation 100)

### 100. Roche Lobe Radius
```
R_L = 0.49q^{2/3} a / [0.6q^{2/3} + ln(1+q^{1/3})]
```
q = M_p/M_*

---

## UQFF Integration Notes

All 100 equations can be integrated into the UQFF framework via:

1. **Buoyancy Modulation (F_U_Bi_i):** Add vacuum density differentials to stabilize/destabilize dynamics
2. **Magnetism Enhancement (Um):** Time-dependent vacuum-modulated magnetic terms
3. **Gravity Duality (Ug):** Primary (superconductivity) vs Newtonian forms
4. **Inertia Operators (Ui):** Gyroscopic torque nullification

### Key UQFF Calibration Constants:
- κ = 0.0005/day
- [SSq] = 0.57
- H_SCm ≈ 0.99
- U_UA ≈ 0.0001
- k_η = 10^{-113}
- β_i ≈ 0.603
- ρ_vac,[UA] = 7.09×10^{-36} J/m³
- ρ_vac,[SCm] = 7.09×10^{-37} J/m³

---

## References

1. Schmidt et al. (2016) - Collision Dynamics of Alpha-Conjugate Nuclei (DOI: 10.1393/ncc/i2016-16394-6)
2. Widom-Larsen (2008/2010) - LENR Theory (arXiv, Pramana)
3. Bond et al. (1991) - Extended Press-Schechter
4. Blandford-Znajek - Jet Power Extraction
5. Sedov-Taylor - Blast Wave Dynamics
6. Tolman-Oppenheimer-Volkoff - Neutron Star Structure

---

*Last Updated: March 1, 2026*
*Extracted from Grok AI Conversation (August 22, 2025)*
