"""S582: Lunar synodic month 29.53 days ; D_c + D_p - F_TRZ*D_phys - F_TRZ*Phi_res + F_TRZ^2*K_Mex"""
from fractions import Fraction
F=Fraction(1,10); Ph=Fraction(5,6); K=Fraction(25,12); D_p=4; D_c=26
v=D_c+D_p-F*D_p-F*Ph+F*F*K
target=29.53
print(f"Lunar month: {float(v):.6f} vs {target} -> {abs(float(v)-target)/target*100:.4f}%")
