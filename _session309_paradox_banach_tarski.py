"""
S309  --  BANACH-TARSKI PARADOX (1924)

In ZFC set theory + axiom of choice, a solid ball in R^3 can be
decomposed into finitely many disjoint pieces and reassembled, using
only rigid motions, into TWO solid balls each identical to the
original.  Volume is not conserved.  Paradox.

UQFF closure: physical space has a discrete BSFG-cell minimum
volume V_min ~ ell_P^3 / D_BSFG^{3/2}.  Non-measurable sets cannot
exist below V_min.  B-T decomposition uses non-measurable pieces;
hence FORBIDDEN in physical space.  B-T is a purely set-theoretic
artifact with no physical realization.
"""
import math

F_TRZ, Phi_res, K_Mex = 0.1, 5.0/6.0, 25.0/12.0
D_phys, D_BSFG, D_crit = 4, 6, 26

ell_P = 1.616255e-35   # Planck length [m]

print("="*72)
print(" S309  --  BANACH-TARSKI PARADOX")
print("="*72)
print()
print(" Theorem (Banach-Tarski 1924, Fund Math 6):")
print("   Let B = {x in R^3 : |x| <= 1}.")
print("   Exists partition B = A_1 u A_2 u ... u A_n  (n = 5 suffices)")
print("   and rigid motions g_1,..,g_n,h_1,..,h_n such that:")
print("     g_1(A_1) u ... u g_n(A_n)  =  B")
print("     h_1(A_1) u ... u h_n(A_n)  =  B")
print("     (B copied twice from itself.)")
print()
print(" Proof requires:  (a) axiom of choice")
print("                  (b) non-measurable Vitali-style sets")
print()
print(" Apparent violation of volume / mass conservation.")
print()
print("-"*72)
print(" UQFF closure: discrete BSFG minimum volume")
print("-"*72)
print()
V_min = (ell_P ** 3) / (D_BSFG ** 1.5)
print(f"   V_min  =  ell_P^3 / D_BSFG^(3/2)")
print(f"          =  ({ell_P:.4e})^3 / 6^1.5")
print(f"          =  {V_min:.4e}  m^3")
print()
print(" Every physical region of space has volume that is an integer")
print(" multiple of V_min (locked BSFG-cell quantization).  Sets of")
print(" volume below V_min are NOT physically realizable.")
print()
print(" Banach-Tarski pieces A_i are NON-MEASURABLE -- they do not")
print(" possess a well-defined volume in the Lebesgue sense, and any")
print(" approximation requires arbitrarily fine ZFC partitions, hence")
print(" sub-V_min subsets.  UQFF forbids these.")
print()
print(" Conclusion: B-T can be PROVEN in pure ZFC, but is")
print(" UNREALIZABLE on the physical space described by UQFF.")
print()

print("-"*72)
print(" Quantitative bound: how 'fine' must a decomposition be to fail?")
print("-"*72)
print()
# Sphere of radius 1 m, volume = 4 pi/3 m^3
V_unit_ball = 4.0 * math.pi / 3.0
N_cells = V_unit_ball / V_min
print(f"   Unit ball:  V  =  {V_unit_ball:.6f}  m^3")
print(f"   Number of BSFG cells in unit ball  =  {N_cells:.3e}")
print()
print(f" To duplicate the ball, B-T would have to split into > {N_cells:.0e}")
print(" cells.  But B-T proof uses only 5 pieces (Robinson 1947 lower")
print(" bound).  These 5 pieces are NECESSARILY non-measurable and")
print(" hence cannot be assembled from BSFG cells.")
print()

print("-"*72)
print(" Solovay model (1970): connection to UQFF")
print("-"*72)
print()
print(" In ZF + Dependent Choice + 'all sets of reals are Lebesgue")
print(" measurable' (Solovay), Banach-Tarski FAILS.")
print()
print(" UQFF is naturally a Solovay-like physical theory: the locked")
print(" BSFG quantization enforces measurability for every physical")
print(" subset of space.  No non-measurable sets exist.  Therefore")
print(" the axiom of choice is replaced (in physics) by Dependent")
print(" Choice, which is sufficient for all of physics-relevant math")
print(" (countable choice, separable Hilbert spaces, etc.) and does")
print(" NOT yield Banach-Tarski.")
print()

print("-"*72)
print(" Falsifiability:")
print("-"*72)
print(" Any reproducible volume-non-conservation in a closed quantum-")
print(" mechanical system would refute UQFF.  None observed.")
print()
print(" Conversely, observation of mass-energy duplication via finite")
print(" rigid-motion partition (e.g. in any classical or quantum")
print(" experiment) would refute the BSFG-cell quantization.  None")
print(" observed.")
print()
print("="*72)
print(" S309 COMPLETE.  Banach-Tarski forbidden by BSFG cell quantization.")
print(f" V_min = ell_P^3 / D_BSFG^(3/2) = {V_min:.3e} m^3.")
print("="*72)
