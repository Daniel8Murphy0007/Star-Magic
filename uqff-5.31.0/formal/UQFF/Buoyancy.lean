/-
  UQFF.Buoyancy — Ubi, Um, and the total assembly F_U.

  Mirrors compute_Ubi_SOURCE4 and compute_Um_SOURCE4 from MAIN_1_CoAnQi.cpp.
  The composite F_U = (Ug1 + Ug2 + Ug3 + Ug4) - Ubi + Um is stated as a
  definition only; no claim is made about its physical realism here.
-/
import Mathlib.Data.Real.Basic
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic
import UQFF.Constants
import UQFF.Ug

namespace UQFF.Buoyancy

open Real UQFF.Constants UQFF.Ug

/-- Buoyancy enhancement of a single Ug component. -/
noncomputable def Ubi
    (Ugi tn Ω_g Mbh dg ε_sw ρ_sw : ℝ) : ℝ :=
  let galactic := Ω_g * (Mbh / dg)
  let enh      := 1 + ε_sw * ρ_sw
  let osc      := Real.cos (Real.pi * tn)
  beta_i * Ugi * galactic * enh * rho_A_SCm * osc

/-- Universal magnetism Um = (M·R²·ω₀) / r³. -/
noncomputable def Um (b : Body) (r : ℝ) : ℝ :=
  (b.M * b.R ^ 2 * b.omega0) / r ^ 3

/-- F_U total assembly. -/
noncomputable def F_U
    (b : Body)
    (r t tn α k1 k2 k3 k4 δdef δsw v_sw R_b P_core
     ρ_v C_conc Ω_g Mbh dg ε_sw ρ_sw : ℝ) : ℝ :=
  let u1  := Ug1 b r t tn α k1 δdef
  let u2  := Ug2 b r t tn k2 δsw v_sw R_b
  let u3  := Ug3 b t k3 P_core
  let u4  := Ug4 t tn ρ_v C_conc α k4
  let sum := u1 + u2 + u3 + u4
  let bi  := Ubi sum tn Ω_g Mbh dg ε_sw ρ_sw
  let um  := Um b r
  sum - bi + um

/-! ### Targeted lemmas (currently `sorry`). -/

theorem Ubi_zero_when_Ugi_zero
    (tn Ω_g Mbh dg ε_sw ρ_sw : ℝ) :
    Ubi 0 tn Ω_g Mbh dg ε_sw ρ_sw = 0 := by
  sorry

theorem Um_pos_of_pos_inputs
    (b : Body) (r : ℝ) (hr : 0 < r) (hω : 0 < b.omega0) :
    0 < Um b r := by
  sorry

end UQFF.Buoyancy
