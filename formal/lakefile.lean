import Lake
open Lake DSL

package «UQFF» where
  -- Formal verification scaffold for UQFF identities and Clay correspondences.
  -- Every claim that has not been proven mechanically is marked `sorry`.
  -- See README.md for the proof-status table.
  leanOptions := #[
    ⟨`pp.unicode.fun, true⟩,
    ⟨`autoImplicit, false⟩
  ]

require mathlib from git
  "https://github.com/leanprover-community/mathlib4.git" @ "v4.11.0"

@[default_target]
lean_lib «UQFF» where
  roots := #[`UQFF]
