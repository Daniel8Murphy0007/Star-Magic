"""S558: Atm pressure 101.325 kPa = SO5^2 + SSq + Phi_res - F_TRZ*Phi_res"""
from fractions import Fraction
F=Fraction(1,10); Ph=Fraction(5,6); SSq=Fraction(57,100); SO=10
v=SO*SO+SSq+Ph-F*Ph
target=101.325
err=abs(float(v)-target)/target*100
print(f"Atm pressure: {float(v):.6f} kPa vs {target} -> {err:.4f}%")
