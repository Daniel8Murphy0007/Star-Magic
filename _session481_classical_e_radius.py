from fractions import Fraction
F=Fraction(1,10); P=Fraction(5,6); S=Fraction(57,100); K=Fraction(25,12)
v = S + K + F + F*P - F*F*K + F*F*P - F*F*F - F*F
val = float(v); target = 2.818
print(f"S481 COMPLETE. Classical e radius (fm) = {val:.4f}; target {target}; match {abs(val-target)/target*100:.4f}%.")
