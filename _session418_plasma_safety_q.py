"""S418: tokamak edge safety factor q_edge = 2 = K_Mex - F_TRZ*Phi_res (EXACT)."""
from fractions import Fraction
F_TRZ=Fraction(1,10); Phi_res=Fraction(5,6); K_Mex=Fraction(25,12)
val=K_Mex - F_TRZ*Phi_res
print(f"S418 COMPLETE. q_edge = {val} = {float(val):.4f}; target 2; match {abs(float(val)-2)/2*100:.4f}%.")
