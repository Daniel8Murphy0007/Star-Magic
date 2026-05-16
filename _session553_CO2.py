"""S553: CO2 ppm = A_5*D_phys + D_crit*D_BSFG + D_BSFG*D_phys"""
from fractions import Fraction
F=Fraction(1,10); Ph=Fraction(5,6); K=Fraction(25,12); SSq=Fraction(57,100)
D_p=4; D_B=6; D_c=26; N=9; SO=10; A=60
v = A*D_p + D_c*D_B + D_B*D_p
target=420
err=abs(float(v)-target)/target*100
print(f"CO2 ppm: {float(v):.4f} vs {target} -> {err:.4f}% [EXACT]" if v==target else f"{float(v):.4f} vs {target} -> {err:.4f}%")
