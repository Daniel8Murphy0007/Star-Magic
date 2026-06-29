/-
  UQFF.Constants — calibrated numerical constants used throughout the framework.

  These values are *empirically calibrated* in the C++ / Python / JavaScript
  implementations (see MAIN_1_CoAnQi.cpp, CondensedPhysics3.py, index.js).
  In this Lean scaffold they are stated as `noncomputable def`s with no claim
  beyond the literal numerical assignment. No physical interpretation is encoded
  at the type level.
-/
import Mathlib.Data.Real.Basic

namespace UQFF.Constants

open Real

/-- SCm vacuum energy density [J/m³]. Calibrated. -/
noncomputable def rho_A_SCm : ℝ := 7.09e-37

/-- UA vacuum energy density [J/m³]. Calibrated. -/
noncomputable def rho_A_UA : ℝ := 7.09e-36

/-- Buoyancy coupling. Empirical fit, dimensionless. -/
noncomputable def beta_i : ℝ := 0.603

/-- Kappa decay rate [per day]. Empirical fit. -/
noncomputable def kappa : ℝ := 0.0005

/-- [SSq] dimensionless coupling. Empirical fit. -/
noncomputable def SSq : ℝ := 0.57

/-- H_SCm dimensionless. Empirical fit. -/
noncomputable def H_SCm : ℝ := 0.99

/-- U_UA dimensionless. Empirical fit. -/
noncomputable def U_UA : ℝ := 1.0e-4

/-- k_eta. Empirical fit. -/
noncomputable def k_eta : ℝ := 1.0e-113

/-- Speed of light [m/s] (CODATA exact since 1983). -/
noncomputable def c_speed : ℝ := 2.99792458e8

/-- Planck constant [J·s] (CODATA exact since 2019). -/
noncomputable def h_planck : ℝ := 6.62607015e-34

/-- Newton's gravitational constant [m³ kg⁻¹ s⁻²]. CODATA 2022. -/
noncomputable def G_grav : ℝ := 6.67430e-11

/-! ### Sanity statements (currently `sorry`).

These are the *bare minimum* facts a scaffold should eventually prove:
positivity of every empirical constant, and that beta_i lies in (0,1).
None of these are the hard work — they are infrastructure. -/

theorem rho_A_SCm_pos : rho_A_SCm > 0 := by
  unfold rho_A_SCm
  norm_num

theorem rho_A_UA_pos : rho_A_UA > 0 := by
  unfold rho_A_UA
  norm_num

theorem beta_i_in_unit : 0 < beta_i ∧ beta_i < 1 := by
  unfold beta_i
  refine ⟨?_, ?_⟩ <;> norm_num

theorem kappa_pos : kappa > 0 := by
  unfold kappa
  norm_num

theorem SSq_in_unit : 0 < SSq ∧ SSq < 1 := by
  unfold SSq
  refine ⟨?_, ?_⟩ <;> norm_num

end UQFF.Constants
