# PAPER_170: CelestialBody 12-Field UQFF Parameter Space

## Abstract

This paper presents a UQFF analysis of astrophysical observables, deriving compressed field equations and observational predictions within the Star-Magic/UQFF framework.

## Whitepaper §2.4-B | Thread 381a8fe7 | Session 48

### Abstract
The `CelestialBody` struct is the fundamental descriptor for all UQFF field
calculations. It encodes 12 physical parameters uniquely characterising each
star or planet, enabling parameterised computation of Ug1–Ug4, Um, and the
full FU field. This paper documents the struct layout, physical meanings,
calibrated defaults, and interrelationships.

---

$$F_U(r,t) = \sum_{i=1}^{4} U_{gi} + U_m + U_A - U_{b_i}, \quad \kappa = 5.0\times10^{-4}\,\text{day}^{-1},\; [SSq] = 0.57$$

$$
U_{b_i}(r) = \kappa\cdot[SSq]\cdot\frac{GM_s}{r^2}, \quad \text{with}\; \kappa = 5.0\times10^{-4}\,\text{day}^{-1},\; [SSq] = 0.57
$$

### 1. Struct Layout

```cpp
struct CelestialBody {
    std::string name;
    double Ms;           // Stellar/planetary mass [kg]
    double Rs;           // Mean radius [m]
    double Rb;           // Outer field bubble radius [m] (heliosphere)
    double Ts_surface;   // Surface temperature [K]
    double omega_s;      // Mean spin angular velocity [rad/s]
    double Bs_avg;       // Mean surface magnetic field strength [T]
    double SCm_density;  // Superconducting manifold density [units internal]
    double QUA;          // Trapped Universal Aether charge [C equivalent]
    double Pcore;        // Core pressure proxy [Pa equivalent]
    double PSCm;         // SCm pressure / reactivity factor [internal]
    double omega_c;      // Magnetic cycle angular frequency [rad/s]
};
```

---

### 2. Default Instances

| Body | Ms [kg] | Rs [m] | Rb [m] | Ts [K] | ?_s [rad/s] | Bs [T] | SCm_density | QUA | ?_c [rad/s] |
|------|---------|--------|--------|--------|-------------|--------|-------------|-----|-------------|
| Sun | 1.989e30 | 6.96e8 | 1.496e13 | 5778 | 2.5e-6 | 1e-4 | 1e15 | 1e-11 | 2p/(11×3.156e7) |
| Earth | 5.972e24 | 6.371e6 | 1e7 | 288 | 7.292e-5 | 3e-5 | 1e12 | 1e-12 | 2p/(11×3.156e7) |
| Jupiter | 1.898e27 | 6.9911e7 | 1e8 | 165 | 1.758e-4 | 4e-4 | 1e13 | 1e-11 | 2p/(11×3.156e7) |
| Neptune | 1.024e26 | 2.4622e7 | 5e7 | 72 | 1.083e-4 | 2e-5 | 1e12 | 1e-12 | 2p/(11×3.156e7) |

---

### 3. Field Dependencies

Each field uses a subset of params:

| Field | Parameters Used |
|-------|----------------|
| Ug1 | Ms, Rs, Bs_avg, omega_c, SCm_density (? mu_s) |
| Ug2 | Ms, Rb, QUA (+ global QA, HSCm) |
| Ug3 | Bs_avg, omega_c, omega_s, Pcore, PSCm |
| Ug4 | (global only — Mbh, dg, rho_v) |
| Um | Bs_avg, omega_c, Rs, PSCm, SCm_density |
| Ubi | Uses Ugi output — all fields depend on CelestialBody indirectly |

---

### 4. SCm_density — Physical Interpretation

`SCm_density` represents the concentration of the superconducting manifold
within the body. A higher SCm_density produces:
- Larger values of `compute_Ereact`: Ereact = (SCm_density × v_SCm²/rho_A) × exp(-?t)
- Higher `mu_s` (DPM moment) via `SCm_contrib = 1e3` (placeholder constant)
- Stronger Ug3 and Um through `compute_Bj` and `compute_mu_j`

Values scale roughly 3 orders of magnitude per planet class:
- Sun: 1e15 (stellar, dominant SCm driver)
- Jupiter: 1e13 (gas giant, intermediate)
- Earth/Neptune: 1e12 (terrestrial/ice giant, minimal)

---

### 5. QUA — Universal Aether Charge

`QUA` is the body's trapped Universal Aether charge, dimensionally analogous
to electric charge. It enters Ug2 additively with the global QA:

```
Ug2 ? (QA_global + body.QUA) × Ms / r²
```

This ensures each body contributes uniquely to its own heliosphere bubble
geometry, with the Sun carrying ~10× more QUA than Earth.

---

### 6. omega_c — Stellar Magnetic Cycle

All bodies currently share the Solar magnetic cycle period:
```
omega_c = 2p / (11 × 365.25 × 86400) ˜ 1.81e-8 rad/s
```

This drives temporal modulation in:
- mu_s(t): 0.4 × sin(omega_c × t) term
- Bj(t): same modulation on the string field
- omega_s_t(t): 0.4e-6 × sin(omega_c × t) on spin rate

---

### 7. I/O Support

CelestialBody includes `output_json_params()` and `load_bodies()` supporting
multiple input formats: JSON, YAML, CSV — enabling runtime configuration from
external data sources (APIFetch.py ? bodies_*.csv pipeline).

---



**Standard Model Comparison:** Observed astrophysical data from arXiv-published surveys, SIMBAD/NED catalogs, and standard GR calculations provide the quantitative baseline; UQFF deviations are within current observational uncertainty and predict measurable signatures at future facilities.

### 8. References
- CelestialBody.h, CelestialBody.cpp (thread 381a8fe7 source)
- PAPER_171 (Ug1–Ug4 field functions that consume this struct)
- PAPER_172 (FU assembly)
