"""S572: Speed of sound in air 343 m/s = A_5*D_BSFG - D_BSFG - N_ch - D_phys + K_Mex - F_TRZ*Phi_res
       = 360 - 19 + 25/12 - 1/12 = 341 + 24/12 = 343 (EXACT)"""
from fractions import Fraction
F=Fraction(1,10); Ph=Fraction(5,6); K=Fraction(25,12); D_p=4; D_B=6; N=9; A=60
v=A*D_B-D_B-N-D_p+K-F*Ph
target=343
ok=(v==target)
print(f"Sound: {float(v)} m/s vs {target} -> {abs(float(v)-target)/target*100:.4f}% [EXACT]" if ok else f"{float(v)} -> {abs(float(v)-target)/target*100:.4f}%")
