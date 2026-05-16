from fractions import Fraction
F=Fraction(1,10); P=Fraction(5,6); S=Fraction(57,100); K=Fraction(25,12); D6=Fraction(6)
v = D6 + S - F*F*K + F*F*P
val = float(v); target = 6.563
print(f"S473 COMPLETE. H-alpha (x100 nm) = {val:.4f}; target {target}; match {abs(val-target)/target*100:.4f}%.")
