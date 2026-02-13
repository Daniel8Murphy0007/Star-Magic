# Wolfram Functions Quick Reference Guide
**Date:** February 13, 2026  
**Status:** Production-ready (31/31 tests passing)

---

## Quick Start Example

```python
from QCalc import UnifiedFieldSolver, ComputeParams, CONSTANTS

# Create solver
solver = UnifiedFieldSolver()

# Define system parameters (automatic detection: magnetar vs SMBH)
params = ComputeParams(
    query_name="My System",
    M=1.4 * CONSTANTS['M_sun'],  # Mass (kg)
    r=20e3,                       # Radius (m)
    B=1e10,                       # Magnetic field (Tesla)
    t=1e8                         # Time (seconds)
)

# Solve - automatically calls Wolfram functions
result = solver.solve(params)

# Access results
for eq in result['long_form_equations']:
    print(f"{eq['name']}: {eq['result']:.3e} {eq['unit']}")
```

---

## System Detection Logic

The integration automatically determines which Wolfram functions to call based on physical parameters:

| System Type | Detection Criteria | Functions Called |
|-------------|-------------------|------------------|
| **Magnetar** | `0.5 < M/M☉ < 3.0` AND `r < 50 km` | 12 magnetar terms (source14) |
| **SMBH** | `M/M☉ > 10^5` | 15 SMBH terms (source15) |
| **Fallback** | `M` and `r` present | Try both, graceful failure |

**Examples:**
- SGR 0501+4516 (M=1.4 M☉, r=20 km) → **Magnetar** → 12 functions
- Sagittarius A* (M=4.3×10^6 M☉, r=12.7 Mkm) → **SMBH** → 15 functions
- NGC 1365 (M=2×10^6 M☉, r=1 kpc) → **SMBH** → 15 functions (M > 10^4 M☉ fallback)

---

## Available Functions

### SOURCE14 - Magnetar Physics (12 functions)

| Function Name | Description | Key Parameters | Output | Notes |
|--------------|-------------|----------------|--------|-------|
| `calculate_base_gravity_hubble_magnetic` | Base gravity with Hubble expansion + magnetic suppression | M, r, B, tau_B, t | g (m/s²) | Time-dependent, magnetic decay |
| `calculate_uqff_unification_time_reversal` | UQFF unified field with time-reversal factor | Ug1-Ug4, B | F_U (m/s²) | f_TRZ = 0.1 constant |
| `calculate_cosmological_constant_acceleration` | Lambda acceleration term | r | a_Λ (m/s²) | Tiny (~10^-36 m/s²) |
| `calculate_em_acceleration_vacuum_corrected` | EM acceleration with vacuum correction | B, v_surf, tau_B, t | a_EM (m/s²) | Magnetic decay included |
| `calculate_gravitational_wave_spin_down` | GW emission from spin-down | M, r, tau_Omega, t | h_GW (m/s²) | Small (~10^-11 m/s²) |
| `calculate_quantum_uncertainty_heisenberg` | Quantum uncertainty contribution | delta_x, delta_p, psi_integral | a_Q (m/s²) | Tiny (~10^-35 m/s²) |
| `calculate_fluid_density_coupling` | Fluid density coupling (Navier-Stokes) | rho, psi_integral | a_fluid (m/s²) | Density-dependent |
| `calculate_oscillatory_wave_superposition` | Wave superposition (standing + traveling) | M, r, B, t, x | a_wave (m/s²) | Oscillatory (±) |
| `calculate_dark_matter_perturbation` | Dark matter perturbation | M_halo, r | a_DM (m/s²) | Small (<10^10 m/s²) |
| `calculate_magnetic_field_decay` | Exponential B-field decay | B, tau_B, t | B(t) (T) | B(t) = B₀ exp(-t/τ_B) |
| `calculate_spin_evolution_angular_velocity` | Spin decay with time | P, tau_Omega, t | Ω(t) (rad/s) | Ω(t) = Ω₀ exp(-t/τ_Ω) |
| `calculate_time_reversal_factor` | Time-reversal factor constant | B | f_TRZ | Always 0.1 (dimensionless) |

### SOURCE15 - SMBH Physics (15 functions)

| Function Name | Description | Key Parameters | Output | Notes |
|--------------|-------------|----------------|--------|-------|
| `calculate_smbh_time_dependent_mass` | M(t) with accretion | M, M_dot, tau_acc, t | M(t) (kg) | M(t) = M₀(1 + Ṁ₀e^(-t/τ_acc)) |
| `calculate_smbh_base_gravity_mass_evolution` | Base gravity with M(t) | M, M_dot, tau_acc, r, t | g (m/s²) | Uses time-dependent mass |
| `calculate_smbh_uqff_unification` | SMBH UQFF unified field | Ug1-Ug4, B | F_U (m/s²) | Same formula as magnetar |
| `calculate_smbh_cosmological_constant` | Lambda acceleration (SMBH scale) | r | a_Λ (m/s²) | Tiny (~10^-30 m/s²) |
| `calculate_smbh_em_acceleration` | EM acceleration (accretion disk) | B, v_surf, tau_B, t | a_EM (m/s²) | Strong (>10^10 m/s²) |
| `calculate_smbh_gravitational_wave` | GW emission with M(t) | M, M_dot, tau_acc, r, tau_Omega, t | h_GW (m/s²) | Tiny (~10^-10 m/s²) |
| `calculate_smbh_quantum_uncertainty` | Quantum uncertainty (SMBH scale) | delta_x, delta_p, psi_integral | a_Q (m/s²) | Tiny (~10^-40 m/s²) |
| `calculate_smbh_fluid_density` | Fluid density with M(t) | M, M_dot, tau_acc, rho, psi_integral, t | a_fluid (m/s²) | Time-dependent |
| `calculate_smbh_oscillatory_wave_orbital` | Orbital oscillations (light-crossing time) | M, r, t, x | a_wave (m/s²) | Frequency ~ c/r |
| `calculate_smbh_dark_matter_precession` | DM with 30° precession | M_halo, r, precession_angle | a_DM (m/s²) | Factor: sin(30°) = 0.5 |
| `calculate_smbh_magnetic_decay_gauss_conversion` | B-field decay (Gauss→Tesla) | B, tau_B, t | B(t) (T) | Converts 10^4 Gauss to 1 T |
| `calculate_smbh_spin_evolution_relativistic` | Relativistic spin (0.3c initial) | r, tau_Omega, t | Ω(t) (rad/s) | Ω₀ = 0.3c/r |
| `calculate_smbh_precession_factor` | Precession factor sin(θ) | precession_angle | factor | sin(30°) = 0.5 |
| `calculate_smbh_accretion_rate` | Accretion rate decay | M_dot, tau_acc, t | Ṁ(t) | Ṁ(t) = Ṁ₀ exp(-t/τ_acc) |
| `calculate_smbh_schwarzschild_radius` | Schwarzschild radius | M | r_s (m) | r_s = 2GM/c² |

---

## Parameter Requirements

### Required Parameters (Must Be Present)
- `M` - Mass (kg) [used for system detection]
- `r` - Radius/distance (m) [used for system detection]

### Optional Parameters (Wolfram-Specific)
- `B` - Magnetic field strength (Tesla)
- `tau_B` - Magnetic decay time (seconds)
- `tau_Omega` - Spin-down time (seconds)
- `tau_acc` - Accretion timescale (seconds)
- `M_dot` - Dimensionless accretion rate (e.g., 0.01 = 1%)
- `M_halo` - Dark matter halo mass (kg)
- `rho` - Fluid density (kg/m³)
- `P` - Rotation period (seconds)
- `v_surf` - Surface velocity (m/s)
- `delta_x` - Position uncertainty (m)
- `delta_p` - Momentum uncertainty (kg·m/s)
- `psi_integral` - Wavefunction integral (dimensionless, normalized to 1.0)
- `precession_angle` - Precession angle (radians, e.g., 30° = 0.524 rad)

### Time Parameters
- `t` - Evolution time (seconds, default: 0.0)

---

## Default Values & Reference Systems

### Magnetar (SGR 0501+4516)
```python
params = ComputeParams(
    query_name="SGR 0501+4516",
    M=1.4 * CONSTANTS['M_sun'],     # 1.4 solar masses
    r=20e3,                          # 20 km radius
    B=1e10,                          # 10^10 Tesla
    tau_B=4000 * 3.156e7,            # 4000 years
    tau_Omega=10000 * 3.156e7,       # 10,000 years
    P=5.0,                           # 5 second rotation period
    rho=1e17,                        # 10^17 kg/m³ density
    v_surf=1e6,                      # 1,000 km/s surface velocity
    delta_x=1e-3,                    # 1 mm position uncertainty
    delta_p=1e-20,                   # 10^-20 kg·m/s momentum uncertainty
    psi_integral=1.0,                # Normalized wavefunction
    M_halo=1e29,                     # Dark matter halo mass
    t=1e8                            # 100 million seconds (~3 years)
)
```

### SMBH (Sagittarius A*)
```python
params = ComputeParams(
    query_name="Sagittarius A*",
    M=4.3e6 * CONSTANTS['M_sun'],    # 4.3 million solar masses
    r=1.27e10,                       # Schwarzschild radius (~12.7 Mkm)
    B=1e4,                           # 10^4 Gauss (1 Tesla)
    tau_B=1e6 * 3.156e7,             # 1 Myr magnetic decay
    tau_Omega=1e9 * 3.156e7,         # 1 Gyr spin-down
    tau_acc=9e9 * 3.156e7,           # 9 Gyr accretion timescale
    M_dot=0.01,                      # 1% dimensionless accretion rate
    rho=1e-10,                       # Low-density accretion disk
    v_surf=1e5,                      # 100 km/s accretion velocity
    delta_x=1e6,                     # 1,000 km uncertainty
    delta_p=1e-15,                   # Momentum uncertainty
    psi_integral=1.0,                # Normalized wavefunction
    M_halo=4.3e4 * CONSTANTS['M_sun'], # 1% DM halo
    precession_angle=30.0 * np.pi / 180,  # 30° in radians
    t=1e12                           # 1 trillion seconds (~31,000 years)
)
```

---

## Typical Output Values

### Magnetar (M=1.4 M☉, r=20 km, B=10^10 T, t=100 Myr)
```
Base Gravity (Hubble + Magnetic): 4.645e+11 m/s²
UQFF Unification (Time-Reversal): 5.929e+11 m/s²
Cosmological Constant Acceleration: 3.296e-36 m/s²
EM Acceleration (Vacuum Corrected): 1.052e+13 m/s²
Gravitational Wave (Spin-Down): 5.075e-11 m/s²
Quantum Uncertainty (Heisenberg): ~10^-35 m/s²
Fluid Density Coupling: ~10^11 m/s²
Oscillatory Wave Superposition: ±10^12 m/s² (oscillatory)
Dark Matter Perturbation: ~10^9 m/s²
Magnetic Field Decay: 1e10 T → 9.99e9 T (after 4 Myr)
Spin Evolution: Ω(t) = 1.26 rad/s → 1.26 rad/s (slow decay)
Time-Reversal Factor: 0.1 (constant)
```

### SMBH (M=4.3×10^6 M☉, r=12.7 Mkm, t=1 Tyr)
```
SMBH Time-Dependent Mass: M(t) ≈ M₀(1 + 0.01e^(-t/τ)) kg
SMBH Base Gravity (M(t) Evolution): ~10^6 m/s²
SMBH UQFF Unification: ~10^6 m/s² (with f_TRZ = 0.1)
SMBH Cosmological Constant: ~10^-30 m/s²
SMBH EM Acceleration: >10^10 m/s² (accretion disk)
SMBH Gravitational Wave (M(t)): ~10^-10 m/s²
SMBH Quantum Uncertainty: ~10^-40 m/s²
SMBH Fluid Density (M(t)): ~10^6 m/s²
SMBH Oscillatory Wave (Orbital): ±10^8 m/s² (light-crossing freq)
SMBH Dark Matter (Precession): ~10^5 m/s² × sin(30°) = ~5×10^4 m/s²
SMBH Magnetic Decay: 1 T → 0.99 T (after 1 Myr)
SMBH Spin Evolution (Relativistic): Ω(t) = 0.3c/r → decay
SMBH Precession Factor: 0.5 (sin(30°))
SMBH Accretion Rate: Ṁ(t) = 0.01e^(-t/9 Gyr)
SMBH Schwarzschild Radius: 1.27×10^10 m
```

---

## Error Handling

### Graceful Failure
If a Wolfram function fails (missing parameters, invalid values), the integration continues with other functions:

```python
result = solver.solve(params)

# Check for warnings
if '_wolfram_warning' in result['solutions']:
    print(f"Warning: {result['solutions']['_wolfram_warning']}")

# Results still available for functions that succeeded
equations = result['long_form_equations']
print(f"Total equations computed: {len(equations)}")
```

### Common Failure Reasons
1. **Missing required parameters** - e.g., `tau_B` not provided for magnetic decay
2. **Invalid parameter values** - e.g., negative mass, zero radius
3. **Division by zero** - e.g., r = 0 in gravity calculation
4. **Numerical overflow** - e.g., extremely large magnetic fields

**Recommendation:** Always provide M, r, and B at minimum for best results.

---

## Testing

Run comprehensive unit tests:
```bash
pytest QCalc_test.py -v
```

**Expected output:**
```
31 passed in 0.76s

TestSource14MagnetarPhysics (12 tests) - All passing
TestSource15SMBHPhysics (15 tests) - All passing
TestQCalcIntegration (4 tests) - All passing
```

---

## Direct Function Access

You can also call Wolfram functions directly without going through QCalc.solve():

```python
from IPData import create_manual_input
from QCalc import CONSTANTS
from QCalc_Wolfram_Extensions import calculate_base_gravity_hubble_magnetic

# Create parameters
params = create_manual_input(
    "My System",  # Name (1st positional arg)
    M=1.4 * CONSTANTS['M_sun'],
    r=20e3,
    B=1e10,
    tau_B=4000 * 3.156e7
)

# Call function directly
result = calculate_base_gravity_hubble_magnetic(params, t=1e8)

print(f"{result.name}: {result.result:.3e} {result.unit}")
print(f"LaTeX: {result.latex}")
print(f"Substituted: {result.substituted}")
print(f"Parameters used: {result.parameters_used}")
print(f"Notes: {result.notes}")
```

---

## FAQ

**Q: How do I know which functions were called?**  
A: Check the equation names in `result['long_form_equations']`. Magnetar equations include "Hubble", "Magnetic", "Time-Reversal". SMBH equations include "SMBH" prefix.

**Q: Why are some Wolfram equations missing from my results?**  
A: Either (1) your system doesn't match magnetar or SMBH criteria, (2) missing required parameters, or (3) function failed and was caught by graceful error handling.

**Q: Can I force all 27 functions to run?**  
A: Yes, provide a wide range of parameters (M, r, B, tau_B, tau_Omega, tau_acc, M_dot, rho, P, v_surf, delta_x, delta_p, psi_integral, M_halo, precession_angle). The integration will try all applicable functions.

**Q: What if I only want magnetar or SMBH functions?**  
A: Set parameters to match detection criteria exactly (magnetar: M=1.4 M☉, r=20 km; SMBH: M>10^5 M☉).

**Q: How do I add more Wolfram functions?**  
A: Follow the pattern in `QCalc_Wolfram_Extensions.py`: Create pure physics function, return `EquationResult`, update `_compute_wolfram_physics_terms()` to call it.

---

## Support & Documentation

- **Full Integration Guide:** `WOLFRAM_INTEGRATION_COMPLETE.md`
- **Architecture Rules:** `CondensedPhysics.py` (MANDATORY ARCHITECTURE RULES header)
- **Test Suite:** `QCalc_test.py` (31 tests with examples)
- **Extraction Source:** `MAIN_1_CoAnQi.cpp` (source14_wolfram.cpp + source15_wolfram.cpp)

**Contact:** Daniel T. Murphy (daniel.murphy00@gmail.com)  
**Framework:** UQFF 99.9% Solvability (Star-Magic)  
**Copyright:** © 2025-2026 Daniel T. Murphy - All Rights Reserved

---

*Last Updated: February 13, 2026 - Integration Phase 1 Complete (31/31 tests passing)*
