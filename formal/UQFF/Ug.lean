/-
  UQFF.Ug — definitions of the four gravity-component functions Ug1..Ug4.

  These mirror the canonical C++ implementations in MAIN_1_CoAnQi.cpp
  (namespace SOURCE4, lines ~25623-26026), specifically:
    compute_Ug1_SOURCE4, compute_Ug2_SOURCE4,
    compute_Ug3_SOURCE4, compute_Ug4_SOURCE4.

  Layer-summed (26-layer) variants from index.js are left for a later commit;
  this file fixes only the single-layer canonical forms.

  No theorem in this file claims a physical fact. The only claims are
  algebraic identities (e.g. positivity of factors, dimensional consistency
  encoded as monotonicity in the relevant variable).
-/
import Mathlib.Data.Real.Basic
import Mathlib.Analysis.SpecialFunctions.Exp
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic
import UQFF.Constants

namespace UQFF.Ug

open Real UQFF.Constants

/-- Celestial body parameters used by Ug1..Ug4. -/
structure Body where
  M       : ℝ   -- mass [kg]
  R       : ℝ   -- radius [m]
  B0      : ℝ   -- surface magnetic field [T]
  omega0  : ℝ   -- angular frequency [rad/s]
  v       : ℝ   -- characteristic velocity [m/s]
  M_pos   : 0 < M
  R_pos   : 0 < R

/-- Body volume V = (4/3) π R³. -/
noncomputable def vol (b : Body) : ℝ := (4 / 3) * Real.pi * b.R ^ 3

/-- Ug1 — magnetic dipole component.
    Ug1 = k₁ · μ_s · (M / r²) · exp(-α t) · cos(π tₙ) · (1 + δ_def)
    where μ_s = ρ_A · V_body. -/
noncomputable def Ug1
    (b : Body) (r t tn α k1 δdef : ℝ) : ℝ :=
  let mu_s   := rho_A_SCm * vol b
  let gradM  := b.M / r ^ 2
  let dec    := Real.exp (-(α * t))
  let osc    := Real.cos (Real.pi * tn)
  let def_f  := 1 + δdef
  k1 * mu_s * gradM * dec * osc * def_f

/-- Ug2 — charge-reactivity / heliosphere component (step-function gated). -/
noncomputable def Ug2
    (b : Body) (r t tn k2 δsw v_sw R_b : ℝ) : ℝ :=
  let Q_SCm    := rho_A_SCm * vol b
  let Q_UA     := rho_A_UA  * vol b
  let E_react  := rho_A_SCm * v_sw ^ 2 / rho_A_UA * Real.exp (-(kappa * t))
  let S_rb     : ℝ := if r > R_b then 1 else 0
  let sw_fac   := 1 + δsw * v_sw
  k2 * (Q_SCm + Q_UA) * b.M / r ^ 2 * S_rb * sw_fac * H_SCm * E_react

/-- Ug3 — magnetic-string rotation component. -/
noncomputable def Ug3
    (b : Body) (t k3 P_core : ℝ) : ℝ :=
  let rot      := Real.cos (b.omega0 * t * Real.pi)
  let E_react  := rho_A_SCm * b.v ^ 2 / rho_A_UA * Real.exp (-(kappa * t))
  k3 * b.B0 * rot * P_core * E_react

/-- Ug4 — vacuum-concentration component. -/
noncomputable def Ug4
    (t tn rho_v C_conc α k4 : ℝ) : ℝ :=
  k4 * rho_v * C_conc * Real.exp (-(α * t)) * Real.cos (Real.pi * tn)

/-! ### Elementary algebraic facts (currently `sorry`).

Reasonable first targets for this scaffold. -/

theorem vol_pos (b : Body) : 0 < vol b := by
  sorry

theorem Ug1_zero_when_k1_zero
    (b : Body) (r t tn α δdef : ℝ) :
    Ug1 b r t tn α 0 δdef = 0 := by
  sorry

theorem Ug2_zero_outside_bubble
    (b : Body) (r t tn k2 δsw v_sw R_b : ℝ)
    (h : r ≤ R_b) :
    Ug2 b r t tn k2 δsw v_sw R_b =
      (let Q_SCm   := rho_A_SCm * vol b
       let Q_UA    := rho_A_UA  * vol b
       let E_react := rho_A_SCm * v_sw ^ 2 / rho_A_UA * Real.exp (-(kappa * t))
       k2 * (Q_SCm + Q_UA) * b.M / r ^ 2 * 0 * (1 + δsw * v_sw) * H_SCm * E_react) := by
  sorry

end UQFF.Ug
