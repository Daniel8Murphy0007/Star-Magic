"""S559: Earth albedo 0.30 = 3*F_TRZ (EXACT)"""
from fractions import Fraction
F=Fraction(1,10)
v=3*F
target=Fraction(30,100)
err=abs(float(v-target))/float(target)*100
print(f"Albedo: {float(v):.4f} vs {float(target)} -> {err:.4f}% [EXACT]" if v==target else f"{float(v)} -> {err:.4f}%")
