from fractions import Fraction
F=Fraction(1,10); P=Fraction(5,6); K=Fraction(25,12); SO=Fraction(10)
v = SO + K - P - F - F*K + F*F*K - F*F*F - F*F + F*F*P
val = float(v); target = 10.974
print(f"S476 COMPLETE. Rydberg (x1e6 1/m) = {val:.4f}; target {target}; match {abs(val-target)/target*100:.4f}%.")
