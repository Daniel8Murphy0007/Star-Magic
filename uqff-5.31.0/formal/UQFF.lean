/-
  UQFF — Lean 4 formal-verification scaffold (root module)

  This module re-exports the four sub-libraries:
    * Constants    — calibrated numerical constants (rho_A, beta_i, kappa, [SSq], H_SCm)
    * Ug           — definitions of Ug1, Ug2, Ug3, Ug4
    * Buoyancy     — Ubi, Um, F_U total assembly
    * Millennium   — Clay Millennium problem statements + UQFF correspondences

  STATUS: every `theorem` in this scaffold is currently `sorry`.
  See README.md for the per-claim proof-status table and the rules used to decide
  what is publishable, what is heuristic, and what is genuinely open.
-/
import UQFF.Constants
import UQFF.Ug
import UQFF.Buoyancy
import UQFF.Millennium
