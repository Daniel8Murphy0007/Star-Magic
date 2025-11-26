# Side-by-Side Comparison: source69_wolfram.cpp vs source70_wolfram.cpp

## Overview
- **source69_wolfram.cpp**: UQFFCompressionModule (Classes 640-649) - Multi-system compression framework
- **source70_wolfram.cpp**: M51UQFFModule (Classes 650-659) - Whirlpool Galaxy M51 specific physics

---

## Class-by-Class Comparison

| Class ID | source69 (Compression) | source70 (M51 Galaxy) |
|----------|------------------------|----------------------|
| **640/650** | **CompressionExpansionFactorTerm** (cosmology)<br/>• H(t,z) Friedmann expansion<br/>• Ωm=0.3, ΩΛ=0.7, H₀=67.15 km/s/Mpc<br/>• Returns: 1 + H(t,z)·t<br/>• Multi-system cosmology | **M51DipoleMagneticTerm** (magnetic)<br/>• BH magnetic dipole: μ_dipole·B<br/>• μ = I·A·ω_spin<br/>• I=1e20 A, A=1e15 m², ω=1e-4 rad/s<br/>• Central black hole specific |
| **641/651** | **CompressionSuperconductiveCorrectionTerm** (magnetic)<br/>• 1 - B/B_crit (shielding)<br/>• B_crit=1e11 T (magnetar scale)<br/>• Dimensionless correction<br/>• General magnetic suppression | **M51SuperconductorEnergyTerm** (magnetic)<br/>• Ug2 = B_super²/(2μ₀)<br/>• B_super = μ₀·H_aether<br/>• H_aether=1e-6 A/m<br/>• Energy density formulation |
| **642/652** | **CompressionEnvironmentalForceTerm** (environment)<br/>• Sum of 12 effects: F_wind, F_erode, F_lensing, F_mag, F_decay, F_coll, F_evo, F_merge, F_sf, F_SN, F_rad, F_BH<br/>• Configurable subterm map<br/>• Multi-system environmental modeling | **M51ExternalTidalTerm** (gravity)<br/>• Ug3' = G·M_NGC5195/d²<br/>• Companion galaxy NGC 5195<br/>• M=1e10 M_sun, d=50 kpc<br/>• Single tidal interaction |
| **643/653** | **CompressionMassEvolutionTerm** (dynamics)<br/>• M(t) = M₀·(1 + M_sf(t))<br/>• M_sf = (SFR·t_yr)/M₀<br/>• Star formation growth over time<br/>• General mass evolution | **M51ReactionEnergyTerm** (nuclear)<br/>• Ug4 = k₄·E_react(t)<br/>• E_react = E₀·exp(-λ·t)<br/>• E₀=1e46 J, λ=0.0005<br/>• Nuclear decay specific |
| **644/654** | **CompressionUg1GravityTerm** (gravity)<br/>• G·M/r² (Newtonian)<br/>• Standard gravity base<br/>• Simple formulation<br/>• Universal gravity | **M51InertialVacuumTerm** (vacuum_energy)<br/>• Ui = λ_I·(ρ_SCm/ρ_UA)·ω_i·cos(π·t_n)·(1+F_RZ)<br/>• ρ_SCm=7.09e-37, ρ_UA=7.09e-36 J/m³<br/>• Vacuum oscillations<br/>• Complex vacuum physics |
| **645/655** | **CompressionUg3ExternalGravityTerm** (gravity)<br/>• G·M_ext/r_ext² (e.g., Sgr A*)<br/>• Configurable external mass<br/>• General external influence<br/>• Any external source | **M51SpiralArmWaveTerm** (wave)<br/>• ψ_spiral = A·exp(-r²/(2σ²))·exp(i(m·φ-ω·t))<br/>• Returns \|ψ\|² density<br/>• m=2 arms, σ=1e3 kpc<br/>• Uses std::complex<double><br/>• Spiral structure specific |
| **646/656** | **CompressionUg4SuperconductiveTerm** (gravity)<br/>• Ug4 = Ug1·f_sc<br/>• Superconductive gravity correction<br/>• Simple scaling factor<br/>• General SC correction | **M51StarFormationForceTerm** (stellar)<br/>• F_SF = k_SF·SFR<br/>• k_SF=1e-10<br/>• Star formation feedback<br/>• Active stellar processes |
| **647/657** | **CompressionQuantumWaveTerm** (quantum)<br/>• Three components:<br/>  - Magnetic: q(v×B)<br/>  - Standing: 2A·cos(kx)·cos(ωt)<br/>  - Quantum: (2π/13.8)·A·Re[exp(i(kx-ωt))]<br/>• Multi-component wavefunction<br/>• General quantum interference | **M51TidalForceTerm** (gravity)<br/>• F_tidal = G·M_companion/d²<br/>• Same physics as class 652<br/>• Used for force calculation<br/>• Companion interaction |
| **648/658** | **CompressionDarkMatterPerturbationTerm** (dark_matter)<br/>• (M_vis+M_DM)·(Δρ/ρ + 3GM/r³)<br/>• 15% visible, 85% DM<br/>• Configurable mass ratio<br/>• General DM framework | **M51DarkMatterCurvatureTerm** (dark_matter)<br/>• (M_vis+M_DM)·(Δρ/ρ + 3GM/r³)<br/>• M_vis=1.2e11 M_sun, M_DM=4e10 M_sun<br/>• M51-specific masses<br/>• Same equation, galaxy-specific |
| **649/659** | **CompressionFluidDynamicsTerm** (fluid)<br/>• ρ_fluid·V·g_base<br/>• ρ=1e-20 kg/m³, V=1e3 m³<br/>• Fluid mass contribution<br/>• General fluid physics | **M51QuantumSpiralIntegralTerm** (quantum)<br/>• (ℏ/√(Δx·Δp))·∫\|ψ\|²dV·(2π/t_Hubble)<br/>• Quantum spiral density integral<br/>• Heisenberg uncertainty<br/>• Spiral-specific quantum |

---

## Category Distribution

### source69_wolfram.cpp (Compression)
1. **Cosmology** (1): Expansion factor
2. **Magnetic** (1): Superconductive correction
3. **Environment** (1): 12-component force sum
4. **Dynamics** (1): Mass evolution
5. **Gravity** (3): Ug1, Ug3', Ug4
6. **Quantum** (1): Multi-component wave
7. **Dark Matter** (1): Perturbation + curvature
8. **Fluid** (1): Fluid dynamics

**Total: 10 classes across 8 categories**

### source70_wolfram.cpp (M51)
1. **Magnetic** (2): BH dipole, superconductor energy
2. **Gravity** (3): External tidal, tidal force, general
3. **Nuclear** (1): Reaction energy decay
4. **Vacuum Energy** (1): Inertial vacuum oscillation
5. **Wave** (1): Spiral arm density waves
6. **Stellar** (1): Star formation feedback
7. **Dark Matter** (1): Curvature with M51 masses
8. **Quantum** (1): Spiral integral

**Total: 10 classes across 8 categories**

---

## Key Differences

### Approach Philosophy
| Aspect | source69 (Compression) | source70 (M51) |
|--------|------------------------|----------------|
| **Scope** | Multi-system framework | Single galaxy system |
| **Applicability** | 19+ astrophysical objects | M51 + NGC 5195 only |
| **Parameters** | Configurable/generic | M51-specific constants |
| **Focus** | General compression physics | Galaxy interaction dynamics |
| **Environmental** | 12 distinct force components | Tidal + star formation only |

### Physics Emphasis
| Category | source69 | source70 |
|----------|----------|----------|
| **Cosmology** | ✓ Friedmann expansion | ✗ Not included |
| **Magnetic** | ✓ Simple B/B_crit | ✓✓ BH dipole + SC energy |
| **Gravity** | ✓✓✓ Ug1, Ug3', Ug4 | ✓✓✓ External + tidal + general |
| **Nuclear** | ✗ Not included | ✓ Reaction decay |
| **Vacuum** | ✗ Not included | ✓ Inertial oscillation |
| **Wave** | ✗ Not included | ✓ Spiral density waves |
| **Stellar** | ✗ Not included | ✓ SF feedback |
| **Quantum** | ✓ Interference | ✓ Spiral integral |
| **Dark Matter** | ✓ Generic | ✓ M51-specific |
| **Fluid** | ✓ Fluid dynamics | ✗ Not included |
| **Environment** | ✓ 12 components | ✗ Simplified |
| **Dynamics** | ✓ Mass evolution | ✗ Not included |

### Mathematical Complexity
| Feature | source69 | source70 |
|---------|----------|----------|
| **Complex numbers** | No | Yes (spiral waves) |
| **Exponential decay** | No | Yes (reaction energy) |
| **Oscillations** | Yes (quantum wave) | Yes (vacuum, spiral) |
| **Multi-component** | Yes (12 environmental) | Yes (spiral phases) |
| **Time evolution** | Yes (H(t,z), M(t)) | Yes (E_react(t), ψ(t)) |

---

## Unique to source69 (Compression)
1. **Friedmann cosmology** - H(t,z) expansion factor
2. **12-component environmental forces** - Comprehensive force modeling
3. **Star formation mass evolution** - M(t) = M₀(1 + M_sf(t))
4. **Fluid dynamics** - ρ_fluid·V·g_base
5. **Three-part quantum wave** - Magnetic + standing + quantum interference
6. **Generic parameter framework** - Configurable for multiple systems

## Unique to source70 (M51)
1. **Black hole magnetic dipole** - μ = I·A·ω_spin
2. **Superconductor energy density** - B²/(2μ₀) formulation
3. **Nuclear reaction decay** - E₀·exp(-λ·t)
4. **Inertial vacuum oscillations** - ρ_SCm/ρ_UA ratios with cos(πt_n)
5. **Spiral arm density waves** - Complex ψ_spiral with |ψ|²
6. **Star formation feedback** - k_SF·SFR force
7. **Companion galaxy tidal** - NGC 5195 specific
8. **Quantum spiral integral** - (ℏ/√(Δx·Δp))·∫|ψ|²dV·(2π/t_Hubble)
9. **M51-specific constants** - Hard-coded galaxy parameters

---

## Shared Physics Concepts

Both files implement:
- **Gravity terms** (different formulations)
- **Dark matter** (same equation, different parameters)
- **Magnetic effects** (different mechanisms)
- **Quantum physics** (different approaches)

---

## Integration Path

### source69 Usage
```cpp
// Apply to any system: Magnetar, Sgr A*, Pillars, etc.
CompressionExpansionFactorTerm c640;
double expansion = c640.compute(t, z, params);

CompressionEnvironmentalForceTerm c642;
c642.setSubterm("F_wind", 1e-10);
c642.setSubterm("F_erode", 5e-11);
double F_env = c642.compute(t, params);
```

### source70 Usage
```cpp
// M51-specific calculations
M51DipoleMagneticTerm c650;
double ug1 = c650.compute(B_M51, params);

M51SpiralArmWaveTerm c655;
double psi_sq = c655.compute(r, phi, t, params); // |ψ|²
```

---

## Wolfram Export Functions

Both files provide:
- `exportSource69ToWolfram()` - Returns Mathematica code for classes 640-649
- `exportSource70ToWolfram()` - Returns Mathematica code for classes 650-659

Each class implements `toWolfram()` method for individual export.

---

## Systems Coverage

### source69 Systems (19+)
- Magnetar SGR 1745-2900
- Sagittarius A*
- Pillars of Creation (M16)
- Westerlund 2
- Tapestry of Blazing Starbirth
- NGC 2525, NGC 3603
- Bubble Nebula
- Antennae Galaxies
- Horsehead Nebula
- NGC 1275, NGC 1792
- Hubble Ultra Deep Field
- *...and more*

### source70 Systems (2)
- Whirlpool Galaxy M51
- NGC 5195 (companion)

---

## Conclusion

**source69** provides a **general-purpose compression framework** applicable to multiple astrophysical systems with configurable parameters and comprehensive environmental modeling.

**source70** delivers **M51-specific physics** with specialized terms for spiral galaxy dynamics, black hole interactions, and tidal effects from a companion galaxy.

**Together**, they represent complementary approaches: **breadth** (source69) vs **depth** (source70) in UQFF physics modeling.

---

**Total Classes: 659**  
**Last Update: November 25, 2025**  
**Copyright: Daniel T. Murphy**
