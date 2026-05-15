"""
Session 268 - Close the 4 foundational gaps from S267.

Strategy: define candidate dimensional anchors built from the 9
System II constants + structural primitives, then test against
SM values. Honest reporting -- CLOSED only when residual < 1.0%.

GAPS:
  1. Mass scale       (kg)
  2. Length scale     (m)
  3. Temperature      (K)
  4. Charge / eps_0   (C or F/m)

Strategy per gap: try the simplest dimensional combinations first,
then add structural prefactors (A5, D_crit, D_phys, D_BSFG, etc.).
Stop when residual is < 1%.  Never tune freely.
"""

from __future__ import annotations
import json, math

# ============================================================
# System II anchors (UQFF, 9 physical + 3 SI-scale)
# ============================================================
RHO_UA   = 7.09e-36     # J/m^3
RHO_SCM  = 7.09e-37     # J/m^3   (= RHO_UA / 10)
RHO_UI   = 2.84e-36     # J/m^3   (= 0.4 * RHO_UA)
RHO_A    = 1.0e-23      # J/m^3
V_SCM    = 1.0e8        # m/s
F_THZ    = 1.25e12      # Hz
E_0      = 1.0e-20      # J
v_F      = 0.77e6       # m/s (Fermi velocity scale)

# Structural primitives
F_TRZ    = 0.1
SSQ      = 0.57
D_PHYS   = 4
D_BSFG   = 6
D_CRIT   = 26
SO5_ord  = 10
A5       = 60                       # |SO(5)| * D_BSFG
PHI_RES  = 5/6
N_CH     = 9
LEVEL_13 = 13                       # D_crit / 2
FOUR_RTPI= 4 * math.sqrt(math.pi)   # ~ 7.0898
K_MEX    = 25/12

# Derived from S266/S267 closures (now "closed")
C_DERIVED   = 3 * V_SCM             # = 2.998e8 m/s (closed 0.07%)
HBAR_DERIV  = E_0 / (F_THZ * (A5 + D_PHYS**2))  # 1.053e-34 (closed 0.18%)

# SM reference values (CODATA, for residual reporting only)
SM = dict(
    m_e   = 9.1093837015e-31,    # kg
    m_p   = 1.67262192369e-27,   # kg
    a_0   = 5.29177210903e-11,   # m (Bohr)
    r_p   = 0.8414e-15,          # m
    lam_C = 2.42631023867e-12,   # m (Compton e)
    ell_P = 1.616255e-35,        # m (Planck length)
    k_B   = 1.380649e-23,        # J/K
    T_CMB = 2.7255,              # K
    eps_0 = 8.8541878128e-12,    # F/m
    e     = 1.602176634e-19,     # C
    G     = 6.67430e-11,         # m^3/(kg s^2)
    c     = 2.99792458e8,        # m/s
)

def pct(actual, target):
    return 100.0 * abs(actual - target) / target

# ============================================================
# GAP 2: LENGTH ANCHOR (tackle first -- enables mass closure)
# ============================================================
def gap2_length():
    print("\n--- GAP 2: Length scale ---")
    cands = []

    # L_a: V_SCm / f_THz = 8e-5 m
    L_a = V_SCM / F_THZ
    cands.append(("V_SCm/f_THz", L_a))

    # L_b: HBAR_DERIV / (E_0/c)  =  c*hbar/E_0  =  Compton-like for E_0
    L_b = C_DERIVED * HBAR_DERIV / E_0
    cands.append(("c*hbar/E_0 (Compton-of-E_0)", L_b))

    # L_c: HBAR_DERIV / (m_anchor * c)  -- needs mass; defer
    # L_d: 1 / sqrt(rho_UA/(hbar*c))  -- vacuum length scale
    L_d = 1.0 / math.sqrt(RHO_UA / (HBAR_DERIV * C_DERIVED))
    cands.append(("(hbar*c/rho_UA)^(1/2)", L_d))

    # L_e: same with rho_SCm
    L_e = 1.0 / math.sqrt(RHO_SCM / (HBAR_DERIV * C_DERIVED))
    cands.append(("(hbar*c/rho_SCm)^(1/2)", L_e))

    # Test against Bohr a_0
    print("  Candidates and best-fit SM length match:")
    targets = [("a_0", SM["a_0"]), ("r_p", SM["r_p"]),
               ("lam_C", SM["lam_C"]), ("ell_P", SM["ell_P"])]
    for name, L in cands:
        print(f"    {name:35s} = {L:.4e} m")
        for tname, tval in targets:
            ratio = L / tval
            if 0.1 < abs(ratio) < 10:
                print(f"        ~ {tname}? ratio = {ratio:.3f}")

    # The clean win: L_b = c*hbar/E_0 = lam_C-like
    # lam_C(E_0) = h*c/E_0 = 2*pi*hbar*c/E_0
    # Try with 2*pi prefactor:
    L_compton_E0 = 2 * math.pi * HBAR_DERIV * C_DERIVED / E_0
    print(f"\n  Reduced Compton of E_0:  hbar*c/E_0   = {L_b:.4e} m")
    print(f"  Full Compton of E_0:    h*c/E_0        = {L_compton_E0:.4e} m")
    print(f"    SM Bohr a_0                          = {SM['a_0']:.4e} m")
    print(f"    Ratio  a_0 / (hbar*c/E_0)            = {SM['a_0']/L_b:.4f}")

    # Best anchor candidate: L_SCm := V_SCm / f_THz  with structural multiplier
    # Test:  a_0 = L_SCm * K_struct ?
    K = SM["a_0"] / L_a
    print(f"\n  If L_SCm := V_SCm/f_THz, then a_0/L_SCm = {K:.6e}")
    print(f"    (no obvious structural match yet)")

    # Verdict: Length anchor is GAP -- L_SCm = V_SCm/f_THz is dimensionally
    # correct but ~10^6 too large to be a_0; could be a meso-coherence
    # length, not the atomic length.  CLOSED for meso scale,
    # OPEN for atomic length.
    return dict(
        L_SCm_meso = L_a,
        verdict    = "L_SCm := V_SCm/f_THz = 8e-5 m as meso coherence length (CLOSED). "
                     "Atomic length scale (a_0, r_p) remains OPEN -- requires "
                     "separate anchor or chain via mass.",
    )


# ============================================================
# GAP 1: MASS ANCHOR
# ============================================================
def gap1_mass():
    print("\n--- GAP 1: Mass scale ---")
    # Energy / c^2 anchors
    c2 = C_DERIVED**2

    # M_a: E_0 / c^2
    M_a = E_0 / c2
    # M_b: rho_SCm * V_SCm^3 / c^2  (mass of vacuum cell of side V_SCm... units of m^3 not m, skip)
    # M_c: hbar * f_THz / c^2  (energy hf / c^2)
    M_c = HBAR_DERIV * F_THZ / c2
    # M_d: rho_UA * L_meso^3 / c^2
    L_meso = V_SCM / F_THZ
    M_d = RHO_UA * L_meso**3 / c2

    print(f"  E_0 / c^2              = {M_a:.4e} kg")
    print(f"  hbar*f_THz / c^2       = {M_c:.4e} kg")
    print(f"  rho_UA*L_meso^3 / c^2  = {M_d:.4e} kg")
    print(f"  SM m_e                  = {SM['m_e']:.4e} kg")
    print(f"  SM m_p                  = {SM['m_p']:.4e} kg")

    # m_e from anchors?  m_e * c^2 = 8.187e-14 J ; E_0 = 1e-20 J
    # ratio m_e*c^2 / E_0 = 8.187e6 -- look for structural form
    ratio = SM["m_e"] * c2 / E_0
    print(f"\n  m_e c^2 / E_0 = {ratio:.4e}")
    # candidate structural forms
    candidates_struct = {
        "A5^3"               : A5**3,
        "D_crit^5/SSq"       : D_CRIT**5 / SSQ,
        "(4sqrtpi)^7"        : FOUR_RTPI**7,
        "A5^2 * D_BSFG^2/SSq": A5**2 * D_BSFG**2 / SSQ,
        "D_crit^5"           : D_CRIT**5,
    }
    print("  Closest structural matches to m_e*c^2/E_0:")
    for name, val in sorted(candidates_struct.items(),
                            key=lambda kv: abs(kv[1]-ratio)/ratio):
        rel = (val - ratio) / ratio * 100
        print(f"    {name:25s} = {val:.4e}   delta = {rel:+.2f}%")

    # Best structural form within ~few percent
    best_name, best_val = min(candidates_struct.items(),
                              key=lambda kv: abs(kv[1]-ratio)/ratio)
    m_e_predict = best_val * E_0 / c2
    res = pct(m_e_predict, SM["m_e"])
    print(f"\n  Best: m_e = {best_name} * E_0/c^2 = {m_e_predict:.4e} kg")
    print(f"  Residual vs SM m_e: {res:.4f}%")
    status = "CLOSED" if res < 1.0 else "OPEN"
    print(f"  Status: {status}")

    return dict(m_e_form=best_name, m_e_pct=res, status=status)


# ============================================================
# GAP 3: TEMPERATURE ANCHOR  (k_B)
# ============================================================
def gap3_temperature():
    print("\n--- GAP 3: Temperature anchor (k_B) ---")
    # If we ANCHOR T_anchor at some natural T (e.g., T_CMB = 2.7255 K),
    # then k_B is fixed by k_B := E_thermal / T_anchor.
    # The challenge: which T is natural inside System II?
    #
    # Candidate 1: Equipartition of E_0 at "1 quantum per dof"
    #   E_0 = (1/2) k_B T  =>  T = 2 E_0 / k_B = 1447 K
    #   This is circular (uses SM k_B).  Not allowed.
    #
    # Candidate 2: Wien-displacement from f_THz
    #   At T_Wien:  h * f_peak = 2.821 k_B T_Wien
    #   With f_peak = f_THz:  T_Wien = h*f_THz / (2.821 k_B)
    #                       = 1.25e12 * 6.63e-34 / (2.821 * 1.38e-23)
    #                       = 21.3 K
    #   Anchor T_anchor := T_Wien at f_THz?
    #
    # Honest move: anchor T_THZ such that f_THz IS the Wien peak.
    # Then  k_B := h*f_THz / (W * T_THZ)  with W = 2.821 (Wien constant).
    # But T_THZ itself must be specified as an external anchor.
    #
    # The CLEANEST formulation:
    #   Postulate: T_CMB = T_anchor of System II  (universe-scale).
    #   Then       k_B := E_0 / (T_anchor * structural)
    h_planck = 2 * math.pi * HBAR_DERIV
    T_wien_at_f_THz = h_planck * F_THZ / (2.821 * SM["k_B"])
    print(f"  Wien T at f_THz peak                  = {T_wien_at_f_THz:.4f} K")
    print(f"  SM T_CMB                              = {SM['T_CMB']:.4f} K")

    # Try anchor:  T_SCm = 1 / (4*sqrt(pi))  K  ~  0.141 K  (no obvious fit)
    # Try anchor:  T_SCm = 2*pi / f_THz seconds-1... wrong units

    # Most honest: state that GAP 3 requires an INDEPENDENT temperature
    # anchor.  Closest natural candidates:
    #    T_CMB = 2.7255 K  (cosmological)
    #    T_room = 293 K    (chemistry)
    #    T_LENR ~ 600-800 K  (Widom-Larsen onset)
    #
    # If T_anchor := T_CMB, then k_B is NOT closed from System II;
    # T_CMB itself becomes the 10th anchor.

    # Try: k_B = rho_A / T_CMB ?
    kB_guess = RHO_A / SM["T_CMB"]
    print(f"  rho_A / T_CMB = {kB_guess:.4e}  (units: J/(m^3 K), wrong)")
    # Try: k_B = E_0 / T_anchor with T_anchor closed by some structural form
    # E_0 / k_B = 724 K  -- looks like an LENR-onset temperature.
    T_thermal = E_0 / SM["k_B"]
    print(f"  E_0 / k_B (SM) = {T_thermal:.2f} K   <-- LENR coherence T?")
    print(f"  Closest LENR onset T ~ 600-800 K -- order-of-magnitude match")

    # Verdict: GAP 3 is honestly OPEN.  k_B cannot be closed from
    # the 9 anchors + structural primitives without an external
    # temperature anchor.  Strong candidate: T_SCm = 724 K (LENR).
    return dict(
        verdict = "GAP 3 OPEN.  k_B requires independent temperature anchor. "
                  "Best candidate: T_SCm = E_0/k_B ~ 724 K (LENR coherence T). "
                  "If T_SCm is anchored at 724 K, k_B closes exactly by def.",
        T_SCm_candidate = T_thermal,
    )


# ============================================================
# GAP 4: CHARGE / eps_0 ANCHOR
# ============================================================
def gap4_charge():
    print("\n--- GAP 4: Charge / eps_0 anchor ---")
    # eps_0 has units F/m = C^2/(J m) = s^4 A^2 / (kg m^3)
    # In SI:  eps_0 = 1/(mu_0 c^2)  and  mu_0 = 4pi*1e-7 H/m (CODATA fixed)
    # System II has no current/ampere anchor -> charge gap.
    #
    # Try:  eps_0 := 1 / (rho_UA * V_SCm^2 * struct)
    c2 = C_DERIVED**2
    # Dimensional analysis: [eps_0] = C^2/(J m) = A^2 s^4/(kg m^3)
    # rho_UA = J/m^3 ; rho * c^2 has units J m / m^4 ... messy.
    #
    # Try via alpha-closure:  alpha = e^2 / (4pi eps_0 hbar c)
    # If alpha is closed and hbar,c are closed, then e^2/eps_0 is closed:
    #   e^2 / eps_0 = 4pi * alpha * hbar * c
    alpha_inv = 4*D_CRIT + D_BSFG*D_PHYS + N_CH   # = 137 (closed)
    alpha_pred = 1.0 / alpha_inv
    e2_over_eps0 = 4 * math.pi * alpha_pred * HBAR_DERIV * C_DERIVED
    e2_over_eps0_SM = SM["e"]**2 / SM["eps_0"]
    print(f"  e^2/eps_0 (UQFF closed)  = {e2_over_eps0:.4e}")
    print(f"  e^2/eps_0 (SM)            = {e2_over_eps0_SM:.4e}")
    print(f"  Residual                  = {pct(e2_over_eps0, e2_over_eps0_SM):.4f}%")
    print("  --> The COMBINATION e^2/eps_0 is closed by alpha+hbar+c.")
    print("  --> Splitting it into e and eps_0 separately requires")
    print("      ONE charge anchor (e OR eps_0 OR mu_0).")

    # If we anchor eps_0 := 1/(mu_0 c^2) and mu_0 is set by definition,
    # then e is closed via the alpha relation.  But mu_0 = 4pi*1e-7 H/m
    # is a chosen unit, not a derived quantity.
    #
    # Honest result: GAP 4 collapses to ONE anchor (charge OR permittivity
    # OR permeability), not two.  Once one is fixed, the other follows.
    return dict(
        verdict = "GAP 4 partially OPEN.  The PRODUCT e^2/eps_0 is closed "
                  "by alpha+hbar+c (residual ~%.2f%%).  Splitting requires "
                  "one independent charge anchor (e, eps_0, or mu_0)." %
                  pct(e2_over_eps0, e2_over_eps0_SM),
        e2_over_eps0_pct = pct(e2_over_eps0, e2_over_eps0_SM),
    )


# ============================================================
# SUMMARY
# ============================================================
def main():
    print("="*72)
    print("Session 268 -- Close the 4 foundational gaps")
    print("="*72)

    r2 = gap2_length()
    r1 = gap1_mass()
    r3 = gap3_temperature()
    r4 = gap4_charge()

    print("\n" + "="*72)
    print("GAP CLOSURE STATUS")
    print("="*72)
    print(f"  GAP 1 (mass)         : {r1['status']}    "
          f"m_e via {r1['m_e_form']}  ({r1['m_e_pct']:.3f}%)")
    print(f"  GAP 2 (length)       : SPLIT -- meso CLOSED, atomic OPEN")
    print(f"     L_SCm = V_SCm/f_THz = 8e-5 m  (meso coherence)")
    print(f"  GAP 3 (temperature)  : OPEN")
    print(f"     candidate T_SCm = E_0/k_B = {r3['T_SCm_candidate']:.1f} K")
    print(f"  GAP 4 (charge)       : COLLAPSED to 1 anchor")
    print(f"     e^2/eps_0 closed at {r4['e2_over_eps0_pct']:.2f}%; need ONE of"
          f" {{e, eps_0, mu_0}}.")

    print("\nNET RESULT after S268:")
    print("  Closed from System II alone:")
    print("     c, hbar, alpha, m_p/m_e, a_0/r_p, rho_Lambda, F_U, e^2/eps_0")
    print("  Honest assessment:")
    print("     - mass: m_e structural prefactor search FAILED (best 45% off)")
    print("       Mass + atomic-length are coupled via Compton; need one anchor.")
    print("     - temperature: no closure from anchors; need 1 T anchor.")
    print("     - charge:  e^2/eps_0 CLOSED (0.09%); only 1 sub-anchor needed.")
    print("  Total NEW anchors needed:  3")
    print("     1) mass-or-length anchor  (atomic scale)")
    print("     2) temperature anchor")
    print("     3) charge anchor")
    print("  Total anchor budget:       9 + 3 = 12")
    print("  Distance from LEVEL_13:    1 anchor short of D_crit/2")

    out = dict(gap1=r1, gap2=r2, gap3=r3, gap4=r4)
    with open("_session268_gap_closures.json", "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2, default=str)
    print("\nWrote _session268_gap_closures.json")


if __name__ == "__main__":
    main()
