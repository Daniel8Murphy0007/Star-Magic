# UQFF Lean 4 Formal-Verification Scaffold

**Status:** Scaffold only. Every non-trivial claim is `sorry`. This directory
exists so future work can attempt mechanical verification of UQFF identities
against the canonical C++ / Python / JavaScript implementations.

This scaffold does **not** prove any Millennium Prize problem, and removing
a `sorry` here is not equivalent to a Clay-acceptable submission. See
`UQFF/Millennium.lean` for the explicit epistemic disclaimers.

## Layout

| File | Purpose |
|---|---|
| `lakefile.lean` | Lake build config, pins Mathlib `v4.11.0`. |
| `lean-toolchain` | Pins Lean toolchain `leanprover/lean4:v4.11.0`. |
| `UQFF.lean` | Root re-export. |
| `UQFF/Constants.lean` | Calibrated constants (`rho_A_SCm`, `beta_i`, `kappa`, `[SSq]`, `H_SCm`). |
| `UQFF/Ug.lean` | `Ug1`..`Ug4` definitions mirroring `compute_Ug{1..4}_SOURCE4` in `MAIN_1_CoAnQi.cpp`. |
| `UQFF/Buoyancy.lean` | `Ubi`, `Um`, `F_U` total assembly. |
| `UQFF/Millennium.lean` | Clay statements + UQFF numerical observations + non-implication stubs. |

## Proof-status table

| Claim | Lean symbol | Status | Honest classification |
|---|---|---|---|
| `rho_A_SCm > 0` | `Constants.rho_A_SCm_pos` | `norm_num` proof | Trivial infrastructure |
| `0 < beta_i < 1` | `Constants.beta_i_in_unit` | `norm_num` proof | Trivial infrastructure |
| `vol b > 0` | `Ug.vol_pos` | `sorry` | Easy, ~5 line proof |
| `Ug1 = 0 when k1 = 0` | `Ug.Ug1_zero_when_k1_zero` | `sorry` | Easy, ~5 line proof |
| `Ug2 step-function gating` | `Ug.Ug2_zero_outside_bubble` | `sorry` | Easy |
| `Ubi = 0 when Ugi = 0` | `Buoyancy.Ubi_zero_when_Ugi_zero` | `sorry` | Easy |
| `Um > 0` for positive inputs | `Buoyancy.Um_pos_of_pos_inputs` | `sorry` | Easy |
| Riemann Hypothesis | `Millennium.RiemannHypothesis` | `sorry` | Open problem (Clay) |
| UQFF Li_26 observation ⇒ RH | `Millennium.uqff_observation_not_imply_RH` | `sorry` | **Likely false as stated** — kept as a permanent guard |
| Yang-Mills mass gap | `Millennium.YangMillsMassGap` | `sorry` (placeholder def) | Open problem (Clay) |
| UQFF F_U ratio ⇒ YM | `Millennium.uqff_observation_not_imply_YM` | `sorry` | Unproven; almost certainly does not hold without a constructive QFT bridge |
| Navier-Stokes regularity | `Millennium.NavierStokesGlobalRegularity` | `True` placeholder | Real formalization needs `Mathlib.Analysis.PDE` |
| Hodge / BSD / P vs NP | placeholders | `sorry` | Need algebraic geometry / number theory / complexity machinery |

## Build (when Lean 4 + elan are installed)

```powershell
cd formal
lake update    # fetches Mathlib v4.11.0 the first time (large download)
lake build
```

The current scaffold is expected to **build with warnings** for every `sorry`.
That is the intended state until specific lemmas are picked up for proof.

## What this scaffold is for

1. **Pin the canonical equation forms.** Lean's type system and `noncomputable def`
   give a tamper-evident record of the exact algebraic shape of `Ug1`..`Ug4`,
   `Ubi`, `Um`, `F_U` matching `MAIN_1_CoAnQi.cpp`.
2. **Separate empirical calibration from logical implication.** Every numerical
   observation (e.g. `Li_26(SSq) ≈ SSq`) lives behind an `axiom`, never a
   `theorem`, so it is impossible to accidentally treat a calibration as a proof.
3. **Hold open a Clay-correspondence boundary.** The `_uqff_observation_not_imply_*`
   theorems sit there as red flags. Closing one requires either (a) an actual
   proof of the Millennium implication or (b) deletion with a recorded reason.

## What this scaffold is NOT

- A submission to Clay Mathematics Institute.
- A claim that UQFF resolves any Millennium problem.
- A peer-reviewed publication. CMI does not accept Lean files; it requires
  publication in a Qualifying Outlet (Annals, Inventiones, JAMS, Acta, Duke,
  CMP, Compositio, etc.) followed by a 2-year wait and community acceptance.

## Recommended next steps

1. Discharge the seven easy infrastructure `sorry`s (`vol_pos`, the four
   "zero-when-zero" lemmas, `Um_pos_of_pos_inputs`).
2. Replace the abstract `riemannZeta` axiom with `Mathlib.NumberTheory.LSeries.RiemannZeta`
   once the build is comfortable importing the full Mathlib zeta module.
3. For Yang-Mills, build a separate `UQFF/Lattice.lean` that formalizes a
   discrete approximation and an explicit continuum-limit estimate — this
   is the **publishable** path identified in the Option 1 survey.
4. Treat each Clay-correspondence `theorem` as a permanent open obligation;
   document any genuine attempt under `formal/notes/`.
