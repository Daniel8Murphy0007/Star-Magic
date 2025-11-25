# source5.cpp - Comprehensive Physics Extraction Analysis

## File Metadata
- **Lines**: 1,063 total
- **Primary Focus**: Modularized UQFF with Self-Expanding Framework 2.0
- **Key Difference from source4**: Separated into classes/modules (CelestialBody, MUGE, FluidSolver, UQFFModule5)
- **Framework**: Enhanced with DarkMatterHaloTerm and VacuumEnergyTerm dynamics

---

## CRITICAL FINDING: SELF-EXPANDING DYNAMICS DISCOVERED

### **UQFFModule5 Class** (lines 150-282) - **CORE DYNAMIC FRAMEWORK**

This is the **self-expanding engine** mentioned in the "look for the dynamics" directive!

#### **Dynamic Capabilities**:
1. **Runtime Term Registration**: `registerDynamicTerm(unique_ptr<PhysicsTerm>)`
2. **Parameter Management**: `setDynamicParameter(name, value)`, `getDynamicParameter(name)`
3. **Dynamic Contributions**: `computeDynamicContributions(params)` - adds to all calculations
4. **State Export**: `exportState(filename)` for cross-module communication
5. **Learning Rate**: `setLearningRate(rate)` for optimization
6. **Enhanced Compute Functions**:
   - `compute_Ug1_enhanced()` = original + dynamic
   - `compute_Ug2_enhanced()` = original + dynamic
   - `compute_Ug3_enhanced()` = original + dynamic
   - `compute_MUGE_enhanced()` = original + dynamic

#### **Key Philosophy** (comments):
```cpp
// All original validated calculations are preserved
// Dynamic terms are additive and optional
```

---

## NEW DYNAMIC PHYSICS TERMS (2 Pre-Built Classes)

### 1. **DarkMatterHaloTerm** (lines 87-117)
```cpp
class DarkMatterHaloTerm : public PhysicsTerm {
private:
    double M_halo;   // Halo mass (kg)
    double r_scale;  // Scale radius (m)
    
public:
    double compute(const map<string, double>& params) const override {
        double r = params.at("r");
        // NFW profile contribution
        double x = r / r_scale;
        double rho_0 = M_halo / (4*PI*r_scale^3*(log(2)-0.5));
        return G * M_halo * log(1+x) / (r*x);
    }
};
```
**Physics**: Navarro-Frenk-White (NFW) dark matter halo density profile
- **Equation**: `g_DM(r) = (G*M_halo*ln(1+r/r_s)) / (r*(r/r_s))`
- **Parameters**: M_halo (halo mass), r_scale (scale radius)
- **Application**: Adds dark matter gravitational contribution at runtime

### 2. **VacuumEnergyTerm** (lines 119-140)
```cpp
class VacuumEnergyTerm : public PhysicsTerm {
private:
    double E_vac_scale;  // Vacuum energy scale
    double lambda;       // Coupling strength
    
public:
    double compute(const map<string, double>& params) const override {
        double t = params.at("t");
        return lambda * E_vac_scale * (1.0 + 0.1*sin(1e-10*t));
    }
};
```
**Physics**: Time-varying vacuum energy fluctuation
- **Equation**: `E_vac(t) = lambda * E_vac_scale * (1 + 0.1*sin(10^-10*t))`
- **Parameters**: E_vac_scale (energy scale), lambda (coupling strength)
- **Application**: Adds oscillating vacuum energy contribution (10% amplitude, ultra-long period)

---

## UQFF MODULE5 ENHANCED METHODS (4 Functions)

### 1. **compute_Ug1_enhanced()** (lines 248-258)
```cpp
double compute_Ug1_enhanced(body, r, t, tn, alpha, delta_def, k1) {
    // Original validated calculation
    double original = compute_Ug1(body, r, t, tn, alpha, delta_def, k1);
    
    // Add dynamic contributions
    map<string, double> params = {{"r", r}, {"t", t}, {"tn", tn}};
    double dynamic = computeDynamicContributions(params);
    
    return original + dynamic;
}
```
**Purpose**: Wraps original Ug1 gravity with dynamic term additions

### 2. **compute_Ug2_enhanced()** (lines 260-271)
Same pattern for Ug2 charge-reactivity gravity

### 3. **compute_Ug3_enhanced()** (lines 273-283)
Same pattern for Ug3 magnetic string rotation

### 4. **compute_MUGE_enhanced()** (lines 285-295)
```cpp
double compute_MUGE_enhanced(sys, res, use_compressed = true) {
    // Original validated calculations
    double original = use_compressed ? compute_compressed_MUGE(sys) 
                                      : compute_resonance_MUGE(sys, res);
    
    // Add dynamic contributions
    map<string, double> params = {{"t", sys.t}, {"r", sys.r}, {"M", sys.M}};
    double dynamic = computeDynamicContributions(params);
    
    return original + dynamic;
}
```
**Purpose**: Enables compressed or resonance MUGE with dynamic term additions

---

## PHYSICS CONSTANTS (Same as source4.cpp)

### Fundamental (3)
1. **PI** = 3.141592653589793
2. **c** = 3.0e8 m/s
3. **G** = 6.67430e-11 m³ kg⁻¹ s⁻²

### Galactic (3)
4. **Omega_g** = 7.3e-16 rad/s
5. **Mbh** = 8.15e36 kg
6. **dg** = 2.55e20 m

### SCm/Aether (6)
7. **v_SCm** = 0.99*c m/s
8. **rho_A** = 1e-23 kg/m³
9. **rho_sw** = 8e-21 kg/m³
10. **v_sw** = 5e5 m/s
11. **QA** = 1e-10 C
12. **Qs** = 0.0

### Decay/Modulation (9)
13-21. **kappa, alpha, gamma, delta_sw, epsilon_sw, delta_def, HSCm, UUA, eta**

### Coupling/Vacuum (7)
22-28. **k1-k4, beta_i, rho_v, C_concentration, f_feedback**

### Navier-Stokes (4)
29. **N** = 32
30. **dt_ns** = 0.1
31. **visc** = 0.0001
32. **force_jet** = 10.0

### Resonance Parameters (15)
33-47. **fDPM, fTHz, Evac_neb, Evac_ISM, Delta_Evac, Fsuper, UA_SCM, omega_i, k4_res, freact, fquantum, fAether, fosc, fTRZ, c_res**

**Total Constants**: 47 (same as source4)

---

## CORE EQUATIONS (Same as source4.cpp)

### Universal Gravity (4 functions)
- `compute_Ug1()` - Dipole-gradient (lines 411-417)
- `compute_Ug2()` - Charge-reactivity (lines 419-425)
- `compute_Ug3()` - Magnetic string rotation (lines 427-432)
- `compute_Ug4()` - Vacuum energy concentration (lines 822-827)

### Universal Buoyancy/Magnetism/Aether (3 functions)
- `compute_Ubi()` - Galactic rotation buoyancy (lines 829-833)
- `compute_Um()` - Billion magnetic strings (lines 434-441)
- `compute_A_mu_nu()` - Metric tensor modulation (lines 835-846)

### Unified Field (1 function)
- `compute_FU()` - Complete unification (lines 848-875)

### Compressed MUGE (9 term functions)
Lines 562-618:
- `compute_compressed_base()` - G*M/r²
- `compute_compressed_expansion()` - 1 + H₀*t
- `compute_compressed_super_adj()` - 1 - B/Bcrit
- `compute_compressed_env()` - Environment factor
- `compute_compressed_Ug_sum()` - SOURCE contributions
- `compute_compressed_cosm()` - Λc²/3
- `compute_compressed_quantum()` - Quantum pressure
- `compute_compressed_fluid()` - Hydrodynamic
- `compute_compressed_perturbation()` - Dark matter + fluctuations
- `compute_compressed_MUGE()` - Full equation (lines 620-637)

### Resonance MUGE (13 term functions + wormhole)
Lines 639-711:
- `compute_aDPM()` - Dipole moment frequency
- `compute_aTHz()` - THz contribution
- `compute_avac_diff()` - Vacuum differential
- `compute_asuper_freq()` - Superconductivity
- `compute_aaether_res()` - Aether resonance
- `compute_Ug4i()` - Reactor gravity
- `compute_aquantum_freq()` - Quantum
- `compute_aAether_freq()` - Aether
- `compute_afluid_freq()` - Fluid
- `compute_Osc_term()` - Oscillation
- `compute_aexp_freq()` - Expansion
- `compute_fTRZ()` - Transition zone
- `compute_a_wormhole()` - Wormhole
- `compute_resonance_MUGE()` - Full equation (lines 713-732)

---

## ASTROPHYSICAL SYSTEMS (7 Systems, same as source4)

Defined in main() (lines 977-984):
1. **SGR 1745-2900** (Magnetar)
2. **Sagittarius A*** (SMBH)
3. **Tapestry of Blazing Starbirth** (Nebula)
4. **Westerlund 2** (Stellar cluster)
5. **Pillars of Creation** (Molecular cloud)
6. **Rings of Relativity** (Cosmological structure)
7. **Student's Guide to the Universe** (Observable universe)

---

## NAVIER-STOKES FLUID SIMULATION (FluidSolver Class)

### Class Definition (lines 754-765)
```cpp
class FluidSolver {
public:
    vector<double> u, v, u_prev, v_prev, dens, dens_prev;
    
    void add_source(), diffuse(), advect(), project(), set_bnd();
    void step(double uqff_g = 0.0);  // UQFF integration here
    void add_jet_force(double force);
    void print_velocity_field();
};
```

### UQFF Integration (lines 849-860)
```cpp
void FluidSolver::step(double uqff_g) {
    // Add UQFF gravity-like force as body force
    for (int i = 1; i <= N; ++i) {
        for (int j = 1; j <= N; ++j) {
            v[IX(i, j)] += dt_ns * uqff_g;
        }
    }
    // Standard NS: diffuse, project, advect, project
}
```

---

## UNIT TESTS (23 Tests, same as source4)

Lines 908-974: All modular MUGE function tests with expected values

---

## HELPER FUNCTIONS (15 Functions)

Lines 309-441:
1. `step_function(r, Rb)` - Heaviside step
2. `compute_Ereact()` - Reactor efficiency with decay
3. `compute_mu_s()` - Magnetic dipole moment
4. `compute_grad_Ms_r()` - Surface gravity
5. `compute_Bj()` - Magnetic string field
6. `compute_omega_s_t()` - Time-varying rotation
7. `compute_mu_j()` - String dipole moment
8. `output_json_params()` - JSON export
9. `load_bodies()` - CSV body loader
10. `load_muge_systems()` - CSV MUGE loader
11. `simulate_quasar_jet()` - NS with UQFF
12. `write_velocity_to_csv()` - Velocity field export
13. `print_summary_stats()` - Min/max/mean statistics

---

## MAIN FUNCTION ENHANCEMENTS

### Command-line Arguments (lines 920-933)
- `--input-bodies <file>`: Load celestial bodies from CSV
- `--input-muge <file>`: Load MUGE systems from CSV
- `--output <file>`: Write velocity field to CSV

### Summary Statistics (new feature)
```cpp
print_summary_stats(fu_values, "FU");
print_summary_stats(compressed_values, "Compressed MUGE");
print_summary_stats(resonance_values, "Resonance MUGE");
```

---

## STRUCTURES

### 1. **CelestialBody** (lines 142-155)
12 parameters: name, Ms, Rs, Rb, Ts_surface, omega_s, Bs_avg, SCm_density, QUA, Pcore, PSCm, omega_c

### 2. **MUGESystem** (lines 38-58)
18 parameters: name, I, A, omega1, omega2, Vsys, vexp, t, z, ffluid, M, r, B, Bcrit, rho_fluid, g_local, M_DM, delta_rho_rho

### 3. **ResonanceParams** (lines 15-33)
15 parameters: fDPM, fTHz, Evac_neb, Evac_ISM, Delta_Evac, Fsuper, UA_SCM, omega_i, k4_res, freact, fquantum, fAether, fosc, fTRZ, c_res

---

## EXTRACTABLE WOLFRAM TERMS

### **High Priority (Dynamic Framework)**:
1. **DarkMatterHaloTerm** - NFW profile dark matter halo
2. **VacuumEnergyTerm** - Time-varying vacuum energy
3. **UQFFModule5EnhancedUg1** - Ug1 with dynamic contributions
4. **UQFFModule5EnhancedUg2** - Ug2 with dynamic contributions
5. **UQFFModule5EnhancedUg3** - Ug3 with dynamic contributions
6. **UQFFModule5EnhancedMUGE** - MUGE with dynamic contributions

### **Same as source4.cpp (19 terms)**:
7-25. UniversalGravity1-4, UniversalBuoyancy, UniversalMagnetism, UniversalAether, UnifiedField, CompressedMUGE, ResonanceMUGE, 7 astrophysical system terms, ReactorEfficiency, NavierStokesQuasarJet

---

## KEY DIFFERENCES FROM source4.cpp

### **Structural**:
1. **Modularization**: Separated into multiple "modules" (CelestialBody.cpp, MUGE.cpp, FluidSolver.cpp, UnitTests.cpp, main.cpp)
   - Note: All compiled into single source5.cpp file (monolithic)
2. **Error Handling**: Added try-catch blocks and runtime_error exceptions
3. **CSV I/O**: Added `write_velocity_to_csv()` function
4. **Statistics**: Added `print_summary_stats()` for analysis

### **Dynamic Framework**:
1. **UQFFModule5 Class**: Complete self-expanding framework
2. **DarkMatterHaloTerm**: NFW halo physics (not in source4)
3. **VacuumEnergyTerm**: Time-varying vacuum energy (not in source4)
4. **Enhanced Functions**: 4 new `_enhanced()` wrappers
5. **State Export**: `exportState()` for module communication

### **Metadata**:
1. **Version Tracking**: `metadata["version"] = "2.0-Enhanced"`
2. **Creation Date**: `metadata["created"] = "2025-11-08"`
3. **Framework Tag**: `metadata["framework"] = "Self-Expanding UQFF"`

---

## WOLFRAM VALIDATION QUERIES

### **New for source5**:
- `WolframAlpha: "NFW profile dark matter halo"`
- `WolframAlpha: "Navarro Frenk White density profile"`
- `WolframAlpha: "vacuum energy fluctuation cosmology"`
- `WolframAlpha: "time varying cosmological constant"`

### **Same as source4**:
- All constant validations (G, c, rho_v, etc.)
- All astrophysical object queries (SGR1745, SagA*, etc.)
- All equation queries (NS, Hubble, dipole moment, etc.)

---

## RECOMMENDED WOLFRAM COMPANION FILE

**Filename**: `source5_wolfram.cpp`

**Classes**: 25 total (6 new + 19 from source4)

**New Classes**:
1. `DarkMatterHaloNFWTerm` - NFW halo contribution
2. `VacuumEnergyFluctuationTerm` - Time-varying vacuum
3. `UQFFModule5Ug1EnhancedTerm` - Ug1 + dynamic
4. `UQFFModule5Ug2EnhancedTerm` - Ug2 + dynamic
5. `UQFFModule5Ug3EnhancedTerm` - Ug3 + dynamic
6. `UQFFModule5MUGEEnhancedTerm` - MUGE + dynamic

**From source4** (19 classes): UniversalGravity1-4, UniversalBuoyancy, UniversalMagnetism, UniversalAether, UnifiedField, CompressedMUGE, ResonanceMUGE, 7 astrophysical systems, ReactorEfficiency, NavierStokesQuasarJet

---

## EXTRACTION SUMMARY

**Total Extractable Entities**: 62
- **Constants**: 47 (same as source4)
- **Equations/Functions**: 34 (30 core + 4 enhanced wrappers)
- **Systems**: 7 (same as source4)
- **Classes**: 6 (UQFFModule5, CelestialBody, MUGESystem, ResonanceParams, FluidSolver, + 2 dynamic PhysicsTerms)
- **Tests**: 23 unit tests
- **New Dynamic Terms**: 2 (DarkMatterHalo, VacuumEnergy)

**Total Lines**: 1,063
**Physics Density**: ~4.4% constants, ~52% functions/classes, ~12% tests, ~20% infrastructure, ~12% I/O

**Key Innovation**: **Self-expanding framework with runtime term registration and enhanced compute methods**

---

## NOTES

- **Same Watermark**: © 2025 Daniel T. Murphy
- **Same References**: All May 2025 MUGE/UQFF attachment PDFs
- **Architecture**: Monolithic file with commented "module boundaries" (not true separation)
- **Enhancement Philosophy**: "All original validated calculations are preserved; dynamic terms are additive and optional"
- **Cross-Module Communication**: `exportState()` saves parameters to file for import by other modules
