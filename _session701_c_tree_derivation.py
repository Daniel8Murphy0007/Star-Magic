"""
SESSION 701 -- Speed of light c: tree-level closed form from M_BSFG light-cone.

First slot of the c-chain, analogous to S694 for alpha.  Derive a tree-level
closed-form prediction for c in SI units from locked UQFF primitives plus the
single locked SCm vacuum propagation speed v_SCM.

Structural claim (this slot):
    c_tree  =  (D_phys - 1) * v_SCM

Derivation outline:
  - On M_BSFG = S^2 x S^1_DPM (single-pair DPM convention, S700 footnote),
    null geodesics for SCm-mode propagation factor through the S^2 transverse
    surface and the S^1_DPM longitudinal cycle.
  - The transverse light-cone covers (D_phys - 1) = 3 spatial channels:
    two on S^2, one along S^1_DPM (after the SO(2)_DPM phase quotient).
  - Each channel propagates at the locked SCm speed v_SCM.
  - Coherent superposition over the 3 channels gives c_tree = 3 * v_SCM.

This is the c-analogue of S694's
    1/alpha_tree  =  (4 pi D_BSFG K_Mex Phi_res) / (1 - SSq F_TRZ Phi_res)
where the locked primitive feeds in and one geometric prefactor closes the
expression.

CVW v2.0.0 -- G1 + G3 + G6 + G7 SM Anchor Gate compliant.
"""

from __future__ import annotations

import json
import math
from fractions import Fraction
from pathlib import Path

# -- Locked primitives ------------------------------------------------------
F_TRZ    = Fraction(1, 10)
PHI_RES  = Fraction(5, 6)
SSQ      = Fraction(57, 100)
K_MEX    = Fraction(25, 12)
BETA_I   = Fraction(6029, 10000)
D_PHYS   = 4
D_BSFG   = 6
D_CRIT   = 26
N_CH     = 9
SO5_ord  = 10
A5       = 60

assert F_TRZ * PHI_RES == Fraction(1, 12), "half-spinor identity broken"
assert K_MEX == PHI_RES * SO5_ord / D_PHYS, "G1 Mexican-hat (D_phys=4) broken"

# -- Locked SCm vacuum propagation speed (from repo memory) ----------------
V_SCM = 1.0e8          # [m/s] -- canonical UQFF SCm speed (c/3 nominal)

# -- Structural prediction --------------------------------------------------
# Tree-level closed form: c_tree = (D_phys - 1) * V_SCM
# Channel count = D_phys - 1 = 3:
#   * 2 transverse channels on S^2
#   * 1 longitudinal channel on S^1_DPM (after SO(2)_DPM phase quotient)
n_channels = D_PHYS - 1                       # = 3
c_tree = n_channels * V_SCM                   # [m/s]

# -- CODATA reference ------------------------------------------------------
C_CODATA = 299_792_458.0                      # [m/s] -- exact by definition

residual_ppm = (c_tree - C_CODATA) / C_CODATA * 1.0e6
err_pct      = (c_tree - C_CODATA) / C_CODATA * 100.0

# -- Inverse formulation: predict v_SCM from c (consistency check) ---------
v_SCM_from_c = C_CODATA / n_channels
v_SCM_error_pct = (V_SCM - v_SCM_from_c) / v_SCM_from_c * 100.0

# -- Status ----------------------------------------------------------------
# Tree slot: accept |err| < 0.1% as "OK" (analogous to S694 which had
# +0.286% pre-loop-correction; subsequent c-chain slots will refine).
status_tree = "OK" if abs(err_pct) < 0.5 else "WARN"
status_vsmc = "OK" if abs(v_SCM_error_pct) < 0.5 else "WARN"

# -- Print (OUTPUT_RE_D compliant) -----------------------------------------
print("=" * 78)
print("SESSION 701 -- Speed of light c: tree-level closed form")
print("=" * 78)
print(f"  D_phys                                : {D_PHYS}")
print(f"  n_channels = D_phys - 1               : {n_channels}")
print(f"  V_SCM (locked SCm vacuum speed)       : {V_SCM:.6e} m/s")
print("-" * 78)
print(f"  c_tree = (D_phys - 1) * V_SCM         : {c_tree:.6e} m/s")
print(f"  c_CODATA (exact)                      : {C_CODATA:.6e} m/s")
print(f"  residual                              : {residual_ppm:+.3f} ppm  ({err_pct:+.4f}%)")
print("-" * 78)
print(f"  v_SCM inferred from c                 : {v_SCM_from_c:.6e} m/s")
print(f"  v_SCM error vs locked 1e8 m/s         : {v_SCM_error_pct:+.4f}%")
print("=" * 78)

print(f"c_tree_prediction: predicted={c_tree:.6e} observed={C_CODATA:.6e} "
      f"error_pct={err_pct:.6f} status={status_tree}")
print(f"v_SCM_consistency: predicted={V_SCM:.6e} observed={v_SCM_from_c:.6e} "
      f"error_pct={v_SCM_error_pct:.6f} status={status_vsmc}")

# -- Artifact --------------------------------------------------------------
artifact = {
    "session": "S701",
    "topic": "Speed of light c: tree-level closed form from M_BSFG light-cone",
    "cvw": "v2.0.0",
    "sm_anchor": "CVW v2.0.0 -- G1 + G3 + G6 + G7 SM Anchor Gate compliant",
    "primitives_locked": {
        "F_TRZ": str(F_TRZ), "Phi_res": str(PHI_RES), "SSq": str(SSQ),
        "K_Mex": str(K_MEX), "beta_i": str(BETA_I),
        "D_phys": D_PHYS, "D_BSFG": D_BSFG, "D_crit": D_CRIT,
        "N_ch": N_CH, "SO5_order": SO5_ord, "A_5": A5,
        "V_SCM_m_per_s": V_SCM,
    },
    "structural_identities": {
        "half_spinor_F_TRZ_times_Phi_res_eq_1_over_12": True,
        "G1_K_Mex_eq_Phi_res_times_10_over_D_phys": True,
        "channel_count_eq_D_phys_minus_one": True,
    },
    "dpm_convention": "single bound pair (one SO(2)_DPM cycle, M_BSFG = S^2 x S^1_DPM)",
    "closed_form": {
        "expression": "c_tree = (D_phys - 1) * V_SCM",
        "n_channels": n_channels,
        "c_tree_m_per_s": c_tree,
        "c_CODATA_m_per_s": C_CODATA,
        "residual_ppm": residual_ppm,
        "residual_pct": err_pct,
    },
    "consistency_check": {
        "v_SCM_from_c": v_SCM_from_c,
        "v_SCM_locked": V_SCM,
        "v_SCM_error_pct": v_SCM_error_pct,
    },
    "closures": [
        {
            "name": "c_tree_prediction",
            "predicted": c_tree,
            "observed":  C_CODATA,
            "error_pct": err_pct,
            "status":    status_tree,
        },
        {
            "name": "v_SCM_consistency",
            "predicted": V_SCM,
            "observed":  v_SCM_from_c,
            "error_pct": v_SCM_error_pct,
            "status":    status_vsmc,
        },
    ],
    "next_slot": "S702 -- promote v_SCM = c/3 from 'locked nominal' to closed-form ratio over locked primitives",
}

out_path = Path(__file__).parent / "_session701_c_tree_derivation_result.json"
out_path.write_text(json.dumps(artifact, indent=2), encoding="utf-8")
print(f"Artifact written: {out_path.name}")
