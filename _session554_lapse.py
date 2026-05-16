"""S554: Lapse rate 6.5 K/km = D_BSFG + SSq - F_TRZ*Phi_res"""
from fractions import Fraction
F=Fraction(1,10); Ph=Fraction(5,6); SSq=Fraction(57,100); D_B=6
v=D_B+SSq-F*Ph
target=6.5
err=abs(float(v)-target)/target*100
print(f"Lapse rate: {float(v):.6f} K/km vs {target} -> {err:.4f}%")
