"""S563: g = 9.81 m/s^2 = N_ch + Phi_res - F_TRZ^2 * K_Mex"""
from fractions import Fraction
F=Fraction(1,10); Ph=Fraction(5,6); K=Fraction(25,12); N=9
v=N+Ph-F*F*K
target=9.81
print(f"g: {float(v):.6f} vs {target} -> {abs(float(v)-target)/target*100:.4f}%")
