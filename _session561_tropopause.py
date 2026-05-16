"""S561: Tropopause height 11 km = SO5 + Phi_res + K_Mex*F_TRZ"""
from fractions import Fraction
F=Fraction(1,10); Ph=Fraction(5,6); K=Fraction(25,12); SO=10
v=SO+Ph+K*F
target=11
err=abs(float(v)-target)/target*100
print(f"Tropopause: {float(v):.6f} km vs {target} -> {err:.4f}%")
