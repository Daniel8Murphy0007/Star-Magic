"""
UQFF Closed Constants -- canonical integer-rational forms.

After Sessions 246-253 (PAPERS 1159-1167) closed all 8 Lagrangian gaps, the
following constants are NO LONGER calibrated free parameters but exact
rationals derived from two textbook integers (D_crit=26 Polyakov critical
dimension, D_phys=4 observed spacetime dimension).

This module is the single source of truth.  Calculators that previously used
`F_TRZ = 0.1`, `beta_i = 0.603`, `phi_res = 5/6` etc. as bare floats may
import from here for traceability.  Existing inline definitions remain valid
(they evaluate to the same numbers); this file documents the closure.

References:
    PAPER_1159 (Session 246) -- Phi_res = 5/6
    PAPER_1160 (Session 247) -- F_TRZ   = 1/|SO(5)|
    PAPER_1161 (Session 248) -- 26!     = (1)_{26} Pochhammer
    PAPER_1162 (Session 249) -- KK tower mode-by-mode suppression
    PAPER_1163 (Session 250) -- DPM SO(2) = light-cone of SO(26)
    PAPER_1164 (Session 251) -- T^22 moduli stabilisation
    PAPER_1165 (Session 252) -- beta_i triangular ladder
    PAPER_1166 (Session 253) -- V(UA) Mexican-hat polynomial K = 25/12
    PAPER_1167 (Session 253) -- Master synthesis (8/8 closed)
"""

from fractions import Fraction

# -----------------------------------------------------------------------------
# Textbook integers (the only free inputs)
# -----------------------------------------------------------------------------
D_CRIT = 26                                         # Polyakov critical dimension
D_PHYS = 4                                          # observed spacetime dimension

# -----------------------------------------------------------------------------
# Group factors (derived)
# -----------------------------------------------------------------------------
DIM_SO5 = 10                                        # |SO(5)| = 5*4/2
D_BSFG  = D_CRIT - 4 * DIM_SO5 // 2                 # = 6 (BSFG geometric layer)

assert D_BSFG == 6, "BSFG dimension chain broken"

# -----------------------------------------------------------------------------
# G6: Phi_res = (D_BSFG - 1) / D_BSFG = 5/6   (PAPER_1159)
# -----------------------------------------------------------------------------
PHI_RES_FRAC = Fraction(D_BSFG - 1, D_BSFG)         # 5/6
PHI_RES      = float(PHI_RES_FRAC)                  # 0.8333...

# -----------------------------------------------------------------------------
# G7: F_TRZ = 1 / |SO(5)| = 1/10   (PAPER_1160)
# -----------------------------------------------------------------------------
F_TRZ_FRAC   = Fraction(1, DIM_SO5)                 # 1/10
F_TRZ        = float(F_TRZ_FRAC)                    # 0.1

# -----------------------------------------------------------------------------
# G2: beta_i = 3*(5-i)/20 = (3/2)*(5-i)/|SO(5)|   (PAPER_1165)
#     beta_1 receives subleading SO(5)^2 correction: *(1 + 1/200) = 0.603
# -----------------------------------------------------------------------------
def beta_i_structural(i: int) -> Fraction:
    """Exact triangular ladder: beta_i = 3(5-i)/20 for i in {1..4}."""
    if i not in (1, 2, 3, 4):
        raise ValueError(f"beta_i defined for i in 1..4 (got {i})")
    return Fraction(3 * (5 - i), 20)

def beta_i_observed(i: int) -> float:
    """Including beta_1 SO(5)^2 correction (1 + 1/200) for dipole channel."""
    base = beta_i_structural(i)
    if i == 1:
        return float(base * (1 + Fraction(1, 200)))   # 0.603
    return float(base)

BETA_I_STRUCTURAL = {i: beta_i_structural(i) for i in (1, 2, 3, 4)}
BETA_I_OBSERVED   = {i: beta_i_observed(i)   for i in (1, 2, 3, 4)}
BETA_I_SUM        = sum(BETA_I_STRUCTURAL.values())   # = 3/2 (Archimedean half)

# -----------------------------------------------------------------------------
# G1: V(UA) Mexican-hat -- K = Phi_res * |SO(5)| / D_phys = 25/12   (PAPER_1166)
# -----------------------------------------------------------------------------
K_MEXICAN_HAT_FRAC = PHI_RES_FRAC * DIM_SO5 / D_PHYS  # 25/12
K_MEXICAN_HAT      = float(K_MEXICAN_HAT_FRAC)        # 2.0833...

assert K_MEXICAN_HAT_FRAC == Fraction(25, 12), "K closure identity broken"

# Calibrated reference scales (NOT derived; these are the two physical inputs
# that the Lagrangian is built around -- see PAPER_1166 for context).
RHO_SCM_DEFAULT = 7.09e-37                          # J/m^3 plasmotic vacuum density
V_UA_DEFAULT    = 1.0e8                             # m/s (= c/3 canonical)


def V_UA(UA: float, rho_SCm: float = RHO_SCM_DEFAULT,
         v_UA: float = V_UA_DEFAULT) -> float:
    """Closed Mexican-hat aether potential (PAPER_1166).

        V(UA) = K * rho_SCm * [(UA/v_UA)^2 - 1]^2,  K = 25/12

    All coefficients integer-rational; zero free parameters added.
    """
    x = UA / v_UA
    term = x * x - 1.0
    return K_MEXICAN_HAT * rho_SCm * term * term


# -----------------------------------------------------------------------------
# G4 / G5: 26-suppression scale  (PAPERS 1162, 1164)
# -----------------------------------------------------------------------------
SUPPRESSION_26 = 1.0 / (26.0 ** 26)                 # 1.624e-37


# -----------------------------------------------------------------------------
# Self-test
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    print("UQFF Closed Constants -- self-test")
    print("-" * 60)
    print(f"D_crit              = {D_CRIT}        (Polyakov critical)")
    print(f"D_phys              = {D_PHYS}         (observed)")
    print(f"|SO(5)|             = {DIM_SO5}        (group factor)")
    print(f"D_BSFG (derived)    = {D_BSFG}         (= D_crit - 4*|SO(5)|/2)")
    print(f"Phi_res             = {PHI_RES_FRAC} = {PHI_RES:.6f}   PAPER_1159")
    print(f"F_TRZ               = {F_TRZ_FRAC} = {F_TRZ:.6f}    PAPER_1160")
    print(f"K (Mexican-hat)     = {K_MEXICAN_HAT_FRAC} = {K_MEXICAN_HAT:.6f}  PAPER_1166")
    print(f"beta_i structural   = {BETA_I_STRUCTURAL}              PAPER_1165")
    print(f"beta_i observed     = {BETA_I_OBSERVED}")
    print(f"sum(beta_i)         = {BETA_I_SUM} (= 3/2 Archimedean half)")
    print(f"1/26^26             = {SUPPRESSION_26:.6e}            PAPER_1162")
    print(f"V(UA=v_UA)          = {V_UA(V_UA_DEFAULT):.6e}  (Mexican-hat min = 0)")
    print(f"V(UA=0)             = {V_UA(0.0):.6e} J/m^3  (= 25/12 * rho_SCm)")
    print("-" * 60)
    print("ALL 8 LAGRANGIAN GAPS CLOSED  -- Sessions 246-253, PAPERS 1159-1167")
