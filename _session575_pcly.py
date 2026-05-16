"""S575: parsec / light-year = 3.2616 ; Phi_res*D_phys - Phi_res*F_TRZ + F_TRZ^2*Phi_res + F_TRZ^3*D_phys"""
from fractions import Fraction
F=Fraction(1,10); Ph=Fraction(5,6); D_p=4
v=Ph*D_p-Ph*F+F*F*Ph+F*F*F*D_p
target=3.2616
print(f"pc/ly: {float(v):.6f} vs {target} -> {abs(float(v)-target)/target*100:.4f}%")
