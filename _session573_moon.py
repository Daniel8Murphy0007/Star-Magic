"""S573: Moon distance / Earth radii = 60.336 ; A_5 + F_TRZ*Phi_res*D_phys"""
from fractions import Fraction
F=Fraction(1,10); Ph=Fraction(5,6); A=60; D_p=4
v=A+F*Ph*D_p
target=60.336
print(f"Moon/RE: {float(v):.6f} vs {target} -> {abs(float(v)-target)/target*100:.4f}%")
