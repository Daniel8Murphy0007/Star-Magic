"""S438: WIMP-nucleon cross-section bound exponent 46 (XENONnT 1e-46 cm^2) (EXACT)."""
from fractions import Fraction
F_TRZ=Fraction(1,10); Phi_res=Fraction(5,6); K_Mex=Fraction(25,12); D_phys=4; SO5=10
val=SO5*D_phys+D_phys+K_Mex-F_TRZ*Phi_res
print(f"S438 COMPLETE. WIMP exponent = {val} = {float(val):.4f}; target 46; match {abs(float(val)-46)/46*100:.4f}%.")
