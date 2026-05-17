"""
_K_Mex_REAL_derivation.py
==========================
HONEST SymPy derivation of K_Mex from the Mexican-hat aether potential.

Claim (PAPER_1166):  K_Mex = Phi_res * |SO(5)| / D_phys = (5/6)*10/4 = 25/12

The potential:
    V(UA) = K * rho_SCm * [ (UA/v_UA)^2 - 1 ]^2
          = a0 + a2*UA^2 + a4*UA^4

This script performs the actual variation, computes minima, masses,
and curvature tests SYMBOLICALLY -- no narration, no numerical fitting.

Result printed at the bottom is the TRUTH about what the Mexican-hat
shape does and does not fix about K.
"""

import sympy as sp

UA, v, rho, K = sp.symbols('UA v rho K', positive=True, real=True)

# ---------------------------------------------------------------
# 1. Write the postulated potential
# ---------------------------------------------------------------
V = K * rho * ((UA/v)**2 - 1)**2
V_expanded = sp.expand(V)
print("=" * 70)
print("STEP 1.  Postulated Mexican-hat potential")
print("=" * 70)
print(f"V(UA) = {V}")
print(f"     = {V_expanded}")
print()

# Identify coefficients a0, a2, a4
poly = sp.Poly(V_expanded, UA)
coeffs = {d: c for c, d in zip(poly.all_coeffs()[::-1], range(poly.degree()+1))}
a0 = coeffs.get(0, 0)
a2 = coeffs.get(2, 0)
a4 = coeffs.get(4, 0)
print(f"a0 = coeff(UA^0) = {sp.simplify(a0)}")
print(f"a2 = coeff(UA^2) = {sp.simplify(a2)}")
print(f"a4 = coeff(UA^4) = {sp.simplify(a4)}")
print()

# ---------------------------------------------------------------
# 2. Stationarity:  dV/dUA = 0
# ---------------------------------------------------------------
dV = sp.diff(V, UA)
extrema = sp.solve(dV, UA)
print("=" * 70)
print("STEP 2.  Stationarity dV/dUA = 0  (Euler-Lagrange for static UA)")
print("=" * 70)
print(f"dV/dUA = {sp.simplify(dV)}")
print(f"Solutions: UA = {extrema}")
print()

# ---------------------------------------------------------------
# 3. Mass-squared at minimum:  m^2 = V''(v)
# ---------------------------------------------------------------
d2V = sp.diff(V, UA, 2)
m2 = sp.simplify(d2V.subs(UA, v))
print("=" * 70)
print("STEP 3.  Mass-squared at the minimum UA = v")
print("=" * 70)
print(f"V''(UA)     = {sp.simplify(d2V)}")
print(f"V''(v)      = m^2 = {m2}")
print()

# ---------------------------------------------------------------
# 4. The three "Mexican-hat consistency checks" in PAPER_1166
# ---------------------------------------------------------------
print("=" * 70)
print("STEP 4.  PAPER_1166 'consistency checks' -- do they fix K?")
print("=" * 70)

# 4a. Normalization a2^2 / (4 a4) = a0
check1_lhs = sp.simplify(a2**2 / (4*a4))
check1_rhs = sp.simplify(a0)
check1_holds_for_any_K = sp.simplify(check1_lhs - check1_rhs)
print(f"  Check (a):  a2^2/(4 a4) = {check1_lhs}")
print(f"              a0          = {check1_rhs}")
print(f"              difference  = {check1_holds_for_any_K}")
if check1_holds_for_any_K == 0:
    print("              VERDICT: holds IDENTICALLY for ANY value of K.")
    print("                       Does NOT constrain K.")
print()

# 4b. V_min = 0 at UA = v
V_at_v = sp.simplify(V.subs(UA, v))
print(f"  Check (b):  V(v) = {V_at_v}")
print( "              VERDICT: V(v) = 0 holds for ANY K (geometric property")
print( "                       of the Mexican-hat shape).  Does NOT fix K.")
print()

# 4c. Curvature ratio |a2|/a4 = 2 v^2
ratio = sp.simplify(sp.Abs(a2) / a4)
print(f"  Check (c):  |a2|/a4 = {ratio}")
print(f"              VERDICT: ratio = 2 v^2 holds for ANY K.")
print( "                       Does NOT fix K.")
print()

# ---------------------------------------------------------------
# 5. What WOULD fix K?
# ---------------------------------------------------------------
print("=" * 70)
print("STEP 5.  What additional condition would actually fix K?")
print("=" * 70)
print("Solve V''(v) = m^2 for K, given an INDEPENDENTLY known m^2:")
K_from_mass = sp.solve(sp.Eq(m2, sp.Symbol('m_sq')), K)[0]
print(f"  K = {K_from_mass}")
print()
print("So K is fixed ONLY IF m_UA^2 (the physical mass of the aether scalar)")
print("is known from outside the potential -- e.g. from a measurement, an RG")
print("matching condition, or a one-loop Coleman-Weinberg calculation.")
print()
print("PAPER_1166 plugs m^2 = (50/3) rho/v^2 INTO this equation and recovers")
print("K = 25/12.  But m^2 = (50/3) rho/v^2 is NOT independently derived in")
print("the paper -- it is obtained by ALREADY ASSUMING K = 25/12.  This is")
print("circular.")
print()

# ---------------------------------------------------------------
# 6. Independent test: does the group-theoretic combination
#    Phi_res * |SO(5)| / D_phys = 25/12 arise from any property
#    of the Mexican-hat shape itself?
# ---------------------------------------------------------------
print("=" * 70)
print("STEP 6.  Does Mexican-hat shape force K = Phi_res*|SO(5)|/D_phys?")
print("=" * 70)
Phi_res = sp.Rational(5, 6)
SO5 = 10
D_phys = 4
K_claimed = sp.Rational(Phi_res * SO5, D_phys)
print(f"  Claimed:  K = Phi_res * |SO(5)| / D_phys = ({Phi_res})*({SO5})/{D_phys}")
print(f"          = {sp.Rational(5,6) * 10 / 4}  =  25/12  =  {sp.Rational(25,12)}")
print()
print("  Question: is there any symbolic property of the polynomial")
print(f"            V(UA) = K rho [(UA/v)^2 - 1]^2 that requires K = 25/12?")
print()
print("  Answer (from SymPy variations above):  NO.")
print("  The Mexican-hat shape fixes:")
print("    - minimum location:  UA_min = +/- v")
print("    - minimum value:     V_min  = 0")
print("    - degree-4 polynomial structure")
print("  The Mexican-hat shape does NOT fix:")
print("    - the prefactor K")
print("    - the scalar mass m_UA^2 = 8 K rho / v^2  (depends on K)")
print()

# ---------------------------------------------------------------
# 7. VERDICT
# ---------------------------------------------------------------
print("=" * 70)
print("VERDICT  (SymPy, not narration)")
print("=" * 70)
print()
print("CLAIM (PAPER_1166):  K_Mex = 25/12 is derived from first principles.")
print()
print("REALITY  (this script):")
print("  The Mexican-hat scalar potential")
print("       V(UA) = K rho_SCm [(UA/v_UA)^2 - 1]^2")
print("  has a one-parameter family of physically distinct realizations")
print("  indexed by K > 0.  No internal consistency condition of this")
print("  potential fixes K.")
print()
print("  All three 'consistency checks' in PAPER_1166 are tautologies of the")
print("  Mexican-hat polynomial shape:")
print("     a2^2/(4 a4) = a0     holds for every K")
print("     V(v) = 0             holds for every K")
print("     |a2|/a4 = 2 v^2      holds for every K")
print()
print("  The claim K = 25/12 is therefore an IDENTIFICATION with a group-")
print("  theoretic combination, NOT a derivation from the action.")
print()
print("  To make K = 25/12 a REAL derivation requires ONE of:")
print("    (i)   an independent measurement of m_UA (none in the repo)")
print("    (ii)  a Coleman-Weinberg one-loop calculation that produces 25/12")
print("          from the matter content circulating in the loop")
print("    (iii) an RG matching condition at a UV scale where K is fixed")
print("          by a higher-symmetry theory")
print("    (iv)  a dimensional-reduction argument from 26D -> 4D where K")
print("          drops out of integrating over the compactified manifold")
print()
print("  None of (i)-(iv) is present in PAPER_1166, PAPER_1167, or any")
print("  Lagrangian script reviewed in this thread.")
print()
print("=" * 70)
