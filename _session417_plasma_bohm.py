"""S417: Bohm diffusion prefactor 1/16 = F_TRZ*Phi_res - F_TRZ^2*K_Mex (EXACT)."""
from fractions import Fraction
F_TRZ=Fraction(1,10); Phi_res=Fraction(5,6); K_Mex=Fraction(25,12)
val=F_TRZ*Phi_res - F_TRZ**2*K_Mex
tgt=Fraction(1,16)
print(f"S417 COMPLETE. Bohm prefactor = {val} = {float(val):.6f}; target 1/16 = {float(tgt):.6f}; match {abs(float(val)-float(tgt))/float(tgt)*100:.4f}%.")
