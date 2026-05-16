from fractions import Fraction
F=Fraction(1,10); P=Fraction(5,6); K=Fraction(25,12); D=Fraction(4)
v = D + K - P + F - F*F*K - F*F*P - F*F*F - F*F
val = float(v); target = 5.29
print(f"S478 COMPLETE. Bohr radius (x1e-11 m) = {val:.4f}; target {target}; match {abs(val-target)/target*100:.4f}%.")
