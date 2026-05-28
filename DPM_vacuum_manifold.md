# DPM Vacuum Manifold — Architecture Reference

**Star-Magic UQFF | Session 207+ | Author: Daniel T. Murphy**
**Status:** Canonical — do not alter section titles without updating all cross-references.

---

## 1. What Is the DPM?

The **Di-Pseudo-Monopole (DPM)** is the foundational force-generating structure of the UQFF
(Unified Quantum Field Framework). It is **not a particle** — it is a grinding pair of simulated
magnetic mono-poles whose rotation differential generates **all downstream physics**, including
mass formation, element synthesis, gravity, and cosmological acceleration.

> **DPM IS THE FOUNDATION. GM/r² IS THE LAST OBSERVABLE PROJECTION. NEVER SWAP.**
> — `dpm_helpers.py` docstring (canonical axiom)

### DPM Poles

| Pole | Identity | Rotation | Character |
|------|----------|----------|-----------|
| **North** | SCm (SuperConductive material) | Clockwise (CW) | Massless, highest reactivity with trapped Aether |
| **South** | UA' (Trapped Universal Aether) | Counter-clockwise (CCW) | Forms when SCm encapsulates free UA |

### Grinding Sequence (Mass Production Pipeline)

```
Step 0: SCm contacts free UA  →  Big Bang contact event
Step 1: SCm encapsulates UA   →  UA'  (trapped Aether)
Step 2: CW × CCW grinding     →  UA'' (first excitation)
Step 3: Progressive densification →  UA''', UA''''
Step 4: Maximum metallicity   →  UA''''' (densest superconductive metal)
```

The DPM is the **union** of the two Python modules described in this document:

```
DPM = scm_vacuum_manifold.py  (BASE)
    + ua_vacuum_manifold.py   (SUPERSTRUCTURE)
```

> Note: `dpm_vacuum_manifold.py` is the canonical consolidated root import path.
> The restored legacy modules `scm_vacuum_manifold.py` and `ua_vacuum_manifold.py`
> are available for compatibility, direct SCm/UA inspection, and delegated numeric helper calls when importable.

---

## 2. Two-File DPM Architecture

```
┌──────────────────────────────────────────────────────────────────────────┐
│                      DPM VACUUM MANIFOLD SYSTEM                          │
│                                                                          │
│  ┌──────────────────────────────┐   ┌──────────────────────────────┐   │
│  │   scm_vacuum_manifold.py     │   │   ua_vacuum_manifold.py      │   │
│  │   BASE LAYER (SCm)           │   │   SUPERSTRUCTURE (UA)        │   │
│  │                              │   │                              │   │
│  │  ρ_vac_SCm = 7.09e-37 kg/m³ │◄──│  imports SCm constants &     │   │
│  │  1.25 THz phonon resonance   │   │  functions                   │   │
│  │  F_U_Bi_i_99 (99-term sum)   │   │                              │   │
│  │  127 compute functions        │   │  Defines UA layers on top:   │   │
│  │  VDS, DVP, BSH, 26D          │   │  UA'  = ρ_vac_SCm            │   │
│  │  LENR (individual expts)     │   │  UA'' = ρ_vac_SCm(1+β cos)  │   │
│  │  Cosmology, string theory    │   │  UA'''= +λω_s term           │   │
│  │  Millennium Prize classes    │   │  UA''''= +Δ_UA4 term         │   │
│  │                              │   │                              │   │
│  │  ~2830 lines, 100% tested   │   │  14 unique compute fns        │   │
│  └──────────────────────────────┘   └──────────────────────────────┘   │
│                         ▲                        │                      │
│                         └──────────────────────── ┘                     │
│                            DPM = SCm ⊗ UA                               │
└──────────────────────────────────────────────────────────────────────────┘
```

### Module Roles

| Module | Role | Canonical Path |
|--------|------|----------------|
| `scm_vacuum_manifold.py` | BASE: primordial SCm vacuum physics | root |
| `ua_vacuum_manifold.py` | SUPERSTRUCTURE: UA layers built on SCm | root + `pdf/` |

---

## 3. UA Layer Hierarchy (Canonical Equations)

The four excited states of the UA vacuum build progressively on the SCm base:

$$
\text{UA}' = \rho_{\text{vac,SCm}}
$$

$$
\text{UA}'' = \rho_{\text{vac,SCm}} \left(1 + \beta_i \cos(\pi t_n)\right)
$$

$$
\text{UA}''' = \rho_{\text{vac,SCm}} \left(1 + \beta_i \cos(\pi t_n) + \lambda_i \omega_s\right)
$$

$$
\text{UA}'''' = \rho_{\text{vac,SCm}} \left(1 + \beta_i \cos(\pi t_n) + \lambda_i \omega_s + \Delta_{UA4}\right)
$$

Where:
- $\beta_i = 0.6$ — buoyancy coupling coefficient  
- $\lambda_i = 1.0$ — coupling constant  
- $\omega_s = 2.5 \times 10^{-6}$ rad/s — stellar angular frequency  
- $\Delta_{UA4} = 0.1$ — third excited layer offset (named constant `DELTA_UA_FOURTH`)  
- $t_n$ — dimensionless negative-time parameter  

### Full DPM Buoyancy Expression

$$
F_{U,Bi,i}^{\text{DPM}} = F_{U,Bi,i,99} \times \left(\text{UA}' + \text{UA}'' + \text{UA}''' + \text{UA}''''\right)
$$

### Calibration Ratio (Exact)

$$
\frac{\rho_{\text{vac,UA}}}{\rho_{\text{vac,SCm}}} = 10
$$

This ratio bridges the microscopic LENR scale ($F_{U,Bi,i}$, outside-to-inside buoyancy)
to the macroscopic cosmological scale ($F_{U,Bi}$, inside-to-outside buoyancy).

---

## 4. Master Integral Form

The long-form outside-to-inside buoyancy integral:

$$
F_{U,Bi,i} = \int_0^\infty \left[ -F_0 + \frac{GM}{r^2} + \rho_{\text{SCm}} U_{UA} \cos(\pi t_n) \right] dr
$$

Discrete 99-term computational form (used in `scm_vacuum_manifold.py`):

$$
F_{U,Bi,i,99} = \sum_{k=1}^{99} \left( -\beta_i \cdot U_{g,k} \cdot \cos(\pi t_n) \cdot \frac{M}{r^2} \right)
$$

UA coupling term:

$$
U_i = \lambda_i \cdot \frac{\rho_{\text{vac,SCm}}}{\rho_{\text{vac,UA}}} \cdot \omega_s \cdot \cos(\pi t_n) \cdot (1 + \Delta_{UA4})
$$

Master expression: `master_99 = simplify(F_U_Bi_i_99 + Ui)`

---

## 5. Calibrated Physical Constants

| Constant | Value | Source Module |
|----------|-------|--------------|
| `RHO_VAC_SCM` | 7.09e-37 kg/m³ | `scm_vacuum_manifold.py` |
| `RHO_VAC_UA` | 7.09e-36 kg/m³ | `scm_vacuum_manifold.py` |
| `THZ_PHONON` | 1.25e12 Hz | `scm_vacuum_manifold.py` |
| `BETA_I` | 0.6 | `scm_vacuum_manifold.py` |
| `LAMBDA_I` | 1.0 | `scm_vacuum_manifold.py` |
| `OMEGA_S` | 2.5e-6 rad/s | `scm_vacuum_manifold.py` |
| `SSQ` | 0.57 (Rational 57/100) | `scm_vacuum_manifold.py` |
| `KAPPA` | 0.0005 /day | `scm_vacuum_manifold.py` |
| `DELTA_UA_FOURTH` | 0.1 | `ua_vacuum_manifold.py` |
| `DPM_DENSITY_RATIO` | 10.0 | `ua_vacuum_manifold.py` |
| `E_PHONON` | 8.284e-22 J | `ua_vacuum_manifold.py` |
| `S26_3` | 1.4531e26 | `ua_vacuum_manifold.py` |
| `PHI_RES` | 0.84 | `ua_vacuum_manifold.py` |

---

## 6. How UA Completes the DPM with SCm

The two modules form a **complementary pair** — neither is sufficient alone:

| Capability | `scm_vacuum_manifold.py` | `ua_vacuum_manifold.py` |
|-----------|--------------------------|-------------------------|
| Primordial vacuum base | ✅ ρ_vac_SCm | ❌ (imports from SCm) |
| 1.25 THz phonon resonance | ✅ Full Gaussian derivation | ❌ (imports) |
| F_U_Bi_i_99 (99-term sum) | ✅ sympy + Monte-Carlo | ❌ (imports) |
| VDS / DVP / BSH series | ✅ 26D convergent sums | ❌ |
| Individual LENR compute fns | ✅ parkhomov, pons_fleischmann, holmlid... | ❌ |
| 26D Ramanujan / π-decoder | ✅ | ❌ |
| Cosmology / string theory | ✅ individual functions | ❌ |
| **UA 4-layer DPM structure** | ❌ | ✅ UA'–UA'''' |
| **F_U_Bi_i_DPM** (product) | ❌ | ✅ F_U_Bi_i_99 × UA_total |
| **Phonon linewidth ODE** | ❌ | ✅ ua_solve_phonon_linewidth() |
| **F_U_Bi calibration proof** | ❌ | ✅ ua_fubi_calibration_proof() |
| **DPM cosmological acceleration** | ❌ | ✅ ua_cosmological_acceleration() |
| **Flat rotation curves (DPM)** | ❌ | ✅ ua_rotation_curve_flat() |
| **Hubble tension modulation** | ❌ | ✅ ua_hubble_tension_modulation() |
| **Dark energy substitute** | ❌ | ✅ ua_dark_energy_substitute() |
| **LENR cross-experiment summary** | ❌ | ✅ ua_lenr_comparison() (7 expts) |
| **UA vs Casimir/QED/SED** | ❌ | ✅ ua_casimir_comparison() |
| **UA in string/M-theory** | ❌ | ✅ ua_string_brane_embedding() |

### Why This Separation Works

1. **SCm** handles the fundamental vacuum substrate and all first-principles derivations  
2. **UA** handles the layered excitation hierarchy and multi-scale bridge from LENR → cosmology  
3. **Neither module hardcodes system-specific data** — all physics is parameterised  
4. Any downstream module (`CondensedPhysics4.py`, `DPMCosmologyModule.py`, etc.)  
   imports from one or both as needed  

---

## 7. Codebase File Inventory — All UA / DPM Related Files

### Core DPM Python Modules

The restored root modules `scm_vacuum_manifold.py` and `ua_vacuum_manifold.py` now exist again for compatibility
and direct inspection, while `dpm_vacuum_manifold.py` remains the canonical consolidated entry point.

| File | Role | Key DPM Content |
|------|------|-----------------|
| [`scm_vacuum_manifold.py`](../scm_vacuum_manifold.py) | BASE layer SCm vacuum | ρ_vac_SCm, F_U_Bi_i_99, 127 compute fns, LENR, cosmology |
| [`ua_vacuum_manifold.py`](../ua_vacuum_manifold.py) | SUPERSTRUCTURE UA layers | 4-layer DPM, F_U_Bi_i_DPM, phonon linewidth ODE, cosmo fns |
| [`pdf/scm_vacuum_manifold.py`](scm_vacuum_manifold.py) | Canonical pdf/ sync copy of SCm | Mirror of root (must stay in sync) |
| [`pdf/ua_vacuum_manifold.py`](ua_vacuum_manifold.py) | Canonical pdf/ sync copy of UA | Mirror of root (must stay in sync) |
| [`dpm_helpers.py`](../dpm_helpers.py) | DPM seed helpers | dpm_ug1_seed(), dpm_ug2_shell() — no-G seed equations |
| [`DPMCosmologyModule.py`](../DPMCosmologyModule.py) | Pre-Big Bang 26-center DPM | 26-center formation, Belly Button resonance, inflation force |
| [`et_scm_vacuum.py`](../et_scm_vacuum.py) | E(t) in SCm vacuum | SCm vacuum density evolution, quintessence comparison |
| [`source10_gpu_dpm_atlas.py`](../source10_gpu_dpm_atlas.py) | GPU DPM spectral atlas | GPU-accelerated DPM resonance profiling |

### DPM in CondensedPhysics Suite

| File | DPM Usage |
|------|-----------|
| [`CondensedPhysics4.py`](../CondensedPhysics4.py) | `f_UA_prime` term in buoyancy, `grind_opp()`, DPM ratio computations |
| [`CondensedPhysics3.py`](../CondensedPhysics3.py) | `rho_UA_prime` parameter in cross-section calculations |
| [`CondensedPhysics_InputData.py`](../CondensedPhysics_InputData.py) | `rho_vac_UA_prime` input parameter for all system queries |
| [`CondensedPhysics_OutputData.py`](../CondensedPhysics_OutputData.py) | `rho_vac_UA_prime` stored in output recall datasets |
| [`CondensedPhysics_Validation.py`](../CondensedPhysics_Validation.py) | DPM parameter validation ranges |

### DPM C++ Modules

| File | Role |
|------|------|
| [`Core/dpm_foundation.h`](../Core/dpm_foundation.h) | C++ DPM foundation header |
| [`Core/Modules/DPMModule.cpp`](../Core/Modules/DPMModule.cpp) | Core DPM C++ module |
| [`modules/subterms/DPMModule.h`](../modules/subterms/DPMModule.h) | Subterm DPM header |
| [`modules/subterms/DPMModule.cpp`](../modules/subterms/DPMModule.cpp) | Subterm DPM implementation |

### DPM Documentation

| File | Content |
|------|---------|
| [`DPM_MATHEMATICS.md`](../DPM_MATHEMATICS.md) | Complete DPM mathematics, grinding sequence, 26D tensor |
| [`DPM_vacuum_manifold.md`](../DPM_vacuum_manifold.md) | **This file** — architecture reference |
| [`COMPLETE_UQFF_EQUATIONS_REFERENCE.md`](../COMPLETE_UQFF_EQUATIONS_REFERENCE.md) | Ug1–Ug4, Ubi, Um canonical equations |

### DPM Whitepapers (selected)

| Paper | Title |
|-------|-------|
| PAPER_046 | DPM Cosmology Dark Photon Manifold |
| PAPER_147 | UQFF F_DPM Vortical Resonance — DPM Driver Aether Coupling |
| PAPER_149 | UQFF SgrA* MUGE F_DPM Dominance Extreme Gravity |
| PAPER_179 | Star Magic 5-Chapter Theory DPM Universal Field Taxonomy |
| PAPER_1049 | Source10 GPU DPM Spectral Atlas |
| PAPER_1131 | SCm Vacuum Manifold Primordial First Principle |
| PAPER_1135 | SCm Vacuum Manifold Hub Reactor Validation |

### Other UA/DPM References in Python Codebase

| File | UA/DPM Content |
|------|---------------|
| [`cpp_extracted_constants.py`](../cpp_extracted_constants.py) | `f_UA_prime = 0.999` C++-extracted constant |
| [`pdf/scm_vacuum_manifold.py`](scm_vacuum_manifold.py) | `rho_UA_prime` in `dpm_species_index_acp()` and `spectral_ladder_26state()` |

---

## 8. DPM Data Flow in the Architecture

```
USER QUERY → source2.cpp (PRINCIPAL GUI)
                    │
                    ▼
             APIFetch.py (55 APIs)  →  bodies_YYYYMMDD.csv
                    │
                    ▼ (dataset passed to calculators)
     ┌──────────────────────────────────────────┐
     │         DPM VACUUM MANIFOLD              │
     │                                          │
     │  scm_vacuum_manifold.py   (BASE)         │
     │  ├─ F_U_Bi_i_99 (buoyancy)              │
     │  ├─ VDS, DVP, BSH (26D series)          │
     │  ├─ LENR compute functions               │
     │  └─ Cosmology / string theory fns        │
     │                  +                       │
     │  ua_vacuum_manifold.py   (SUPERSTRUCTURE)│
     │  ├─ UA' → UA'''' (4 excited layers)     │
     │  ├─ F_U_Bi_i_DPM (full product)         │
     │  ├─ Phonon linewidth ODE                 │
     │  └─ Cosmological / multi-expt fns        │
     └──────────────────────────────────────────┘
                    │
                    ▼
          OPData.py → uqff_results.json
                    │
                    ▼
     CondensedPhysics_OutputData.py (RECALL)
```

---

## 9. Usage Pattern

### Import the full DPM system

```python
# BASE layer — all SCm physics
from scm_vacuum_manifold import (
    SSQ, KAPPA, RHO_VAC_SCM, RHO_VAC_UA, THZ_PHONON,
    BETA_I, LAMBDA_I, OMEGA_S,
    F_U_Bi_i_99, monte_carlo_fubi_i,
    compute_F_U_Bi_i_numerical, vds_numerical,
    parkhomov_excess_heat, pons_fleischmann_excess_heat,
)

# SUPERSTRUCTURE — UA layers and DPM bridge functions
from ua_vacuum_manifold import (
    UA_prime, UA_double_prime, UA_triple_prime, UA_quad_prime,
    UA_total, F_U_Bi_i_DPM,
    ua_layer_density, ua_dpm_total_density, ua_dpm_buoyancy_factor,
    ua_calibration_ratio,
    ua_solve_phonon_linewidth, ua_linewidth_convergence,
    ua_lenr_comparison, ua_casimir_comparison, ua_string_brane_embedding,
    ua_cosmological_acceleration, ua_rotation_curve_flat,
    ua_hubble_tension_modulation, ua_dark_energy_substitute,
    ua_fubi_calibration_proof,
    DELTA_UA_FOURTH, DPM_DENSITY_RATIO,
)
```

### Compute the full DPM buoyancy at t_n = 0.25

```python
t_n_val = 0.25
total_ua_density = ua_dpm_total_density(t_n_val)     # sum of 4 UA layers
mean_fubi, _, _ = monte_carlo_fubi_i()                # F_U_Bi_i_99 numerical
F_dpm = mean_fubi * total_ua_density                  # full DPM buoyancy
print(f"F_U_Bi_i_DPM = {F_dpm:.4e}")
```

### Access the DPM calibration proof

```python
proof = ua_fubi_calibration_proof()
# {'rho_vac_SCm': 7.09e-37, 'rho_vac_UA': 7.09e-36,
#  'ratio_UA_over_SCm': 10.0, 'F_U_Bi_i_MC_mean_N': ...,
#  'F_U_Bi_cosmological': ..., 'scale_interpretation': '...'}
```

---

## 10. Maintenance Rules

1. **Always keep `pdf/` in sync** — after any edit to root `ua_vacuum_manifold.py` or `scm_vacuum_manifold.py`, run:
   ```powershell
   Copy-Item ua_vacuum_manifold.py pdf\ua_vacuum_manifold.py -Force
   Copy-Item scm_vacuum_manifold.py pdf\scm_vacuum_manifold.py -Force
   ```

2. **No hardcoded system data** — both modules are pure physics calculators. System-specific data lives in `bodies_*.csv` and `CondensedPhysics_InputData.py`.

3. **No repeated initialization blocks** — `ua_vacuum_manifold.py` has exactly one import block, one constants block, and one `if __name__ == "__main__":` block.

4. **DELTA_UA_FOURTH is a named constant** — never write `0.1` in an equation; always use the constant.

5. **DPM density ratio is DPM_DENSITY_RATIO = 10** — never hardcode `10` in calibration code.

6. **Git commit both .py and pdf/ copy together** — use `--no-verify` on all commits.

---

*Last updated: Session 207+ | Daniel T. Murphy*
