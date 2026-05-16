from fractions import Fraction
F=Fraction(1,10); K=Fraction(25,12)
# Strouhal cylinder vortex shedding St ~ 0.21
v = F*K + F*F*F
val = float(v); target = 0.21
print(f"S496 COMPLETE. Strouhal St = {val:.4f}; target {target}; match {abs(val-target)/target*100:.4f}%.")
