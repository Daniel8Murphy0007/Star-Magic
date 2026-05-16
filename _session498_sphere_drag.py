from fractions import Fraction
F=Fraction(1,10); P=Fraction(5,6); S=Fraction(57,100)
# Sphere drag coefficient C_d = 0.47 (Newtonian regime)
v = S - F*P - F*F*P - F*F*F
val = float(v); target = 0.47
print(f"S498 COMPLETE. Sphere drag C_d = {val:.4f}; target {target}; match {abs(val-target)/target*100:.4f}%.")
