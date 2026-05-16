"""S555: Atmospheric scale height 8.5 km = 2*D_phys + SSq - F_TRZ^2"""
from fractions import Fraction
F=Fraction(1,10); SSq=Fraction(57,100); D_p=4
v=2*D_p+SSq-F*F
target=8.5
err=abs(float(v)-target)/target*100
print(f"Scale height: {float(v):.6f} km vs {target} -> {err:.4f}%")
