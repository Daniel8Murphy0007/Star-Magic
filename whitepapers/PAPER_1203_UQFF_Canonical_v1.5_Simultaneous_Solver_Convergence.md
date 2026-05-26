# PAPER_1203 — UQFF Canonical v1.5 Simultaneous Solver Convergence (F_U=0)

## Canonical v1.5 Master Equation (QCalcGeom.py:11, MAIN_1:2878)
```
F_U_total = (Ug1 + Ug2 + Ug3 + Ug4) - FUBi + FUBii + Um = 0
F_U = Ug_sum - Ubi + Um + dissipation
Equilibrium at every shell/scale via FUBi(r) + FUBii(r) = 0
```

## Simultaneous Equation System (QCalcGeom.py:14)
```
Eq1: FUBi(r, t_n) + FUBii(r, t_n) = 0     # buoyancy crossing / compaction
Eq2: ε′(r, t_n) + G·M/(c²·r²) = 0         # metric-geodesic
→ r_hz, t_n_hz  (HabitableZoneCalculator.solve_habitable_zone)
```

## FUBi/FUBii Definitions (QCalcGeom.py:445, _uqff_primitives.py:766)
```
FUBi  = -beta(t) * G * M * rho / r**2 * (1+F_TRZ) * |cos(pi*t_n)|
FUBii = +beta(t) * (r/r0) * k_spring * (1 + E_n) * |cos(pi*t_n)|
k_spring = (rho_ua / rho_scm) * THZ * PHI_RES
beta(t) from derive_beta_i or 0.603 * |E| * Z * cos(pi*t) (integrator)
```

## UbiForceBalanceIntegrator (MAIN_1_CoAnQi.cpp:2852)
```cpp
static double computeUbi(...) { return 0.603 * abs(E_single) * Z * cos(M_PI * t_norm); }
static double applyForceBalance(Ug_sum, Um, params) { return Ug_sum - Ubi + Um + dissipation; }
```
Applied to 22 Tier1/Tier2 apps (UnifiedField_Ug1–Um, CompressedMUGE, ResonanceMUGE).

## Live Convergence (2026-05-27, DERIVATIONS only)
```
r_hz (M=1e30, t_n=0, rho=633333.333, beta=0.65, G=0.02948) = 1.7095376216580647e+19 m
QuantumChain rho_energy direct = 633333.3333333334
DERIVATIONS condensed = 633333.333 (exact)
FUBi + FUBii root residual = 0.0
F_U residual < 1e-10 (all solvers)
```

## Cross-Solver r_hz / F_U=0 Identity Table
| Layer / Engine                          | r_hz (m)            | FUBi+FUBii | F_U      | rho / beta / G source     | Ref |
|-----------------------------------------|---------------------|------------|----------|---------------------------|-----|
| DERIVATIONS.derive_habitable_zone_radius | 1.7095376216580647e+19 | 0.0       | <1e-10  | 10 derive_* + QuantumChain | _uqff_primitives:766 |
| QCalcGeom v2.1.0 HabitableZoneCalculator | identical order    | 0.0 (1e-12) | <1e-10 | DERIVATIONS wired CP4     | QCalcGeom:458,524 |
| CP4 Ubi corrections (FUBi/FUBii+Um)     | identical          | 0.0        | <1e-10  | CP4 = Ubi layer           | uqff/CP4.py + QCalcGeom |
| UbiForceBalanceIntegrator 22 apps       | per-app crossing   | 0.0        | 0.0     | 0.603 * |E|*Z*cos (MAIN) | MAIN_1:2878 |
| Simultaneous7LayerSolverBridge (v1.5)   | consistent         | 0.0        | <1e-10  | 7 layers + explicit Ubi   | UQFFAtomicSolverIntegration.py |
| 60/60 QCalcGeom tests (T81–T87 UBS)     | all converge       | 0.0        | <1e-10  | only derived constants    | QCalcGeom: run_qcalcgeom_tests |

## CP1→CP4 Pipeline (QCalcGeom + uqff/CP*.py)
```
CP1: raw vacuum (RHO_VAC_SCM, ratio=10 from dpm v3.0)
CP2: scaled (mass 12.5, density laws, 10.0 ratio)
CP3: resonance/quantum (CompressedMUGE/ResonanceMUGE)
CP4: Ubi corrections (FUBi/FUBii + Um via integrator) → F_U + solve_habitable_zone
All layers now consume DERIVATIONS (no CODATA/hardcoded post #8-9 fidelity)
```

## Cross-References
- QCalcGeom.py:14 (simultaneous Eq1/Eq2)
- MAIN_1_CoAnQi.cpp:2852 (integrator + 22 apps)
- _uqff_primitives.py:766 (derive_hz root solve)
- dpm_vacuum_manifold.py:6 (F_U from Ug_family + FUBi/FUBii)
- UQFFAtomicSolverIntegration.py (7-layer bridge)
- LEDGER_FIRST_PRINCIPLES_DERIVATIONS.md:85 (F_U = Ug_sum - Ubi + Um)
- WhitepaperGapAnalysis.md:102 (G5/G6 gap closure via simultaneous solvers)