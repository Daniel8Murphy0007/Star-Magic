"""S560: Greenhouse effect 33 K = D_crit + N_ch - K_Mex + F_TRZ*Phi_res"""
from fractions import Fraction
F=Fraction(1,10); Ph=Fraction(5,6); K=Fraction(25,12); D_c=26; N=9
v=D_c+N-K+F*Ph
target=33
err=abs(float(v)-target)/target*100
print(f"Greenhouse: {float(v):.6f} K vs {target} -> {err:.4f}%")
