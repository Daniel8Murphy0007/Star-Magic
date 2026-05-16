from fractions import Fraction
F=Fraction(1,10); P=Fraction(5,6); K=Fraction(25,12); D4=Fraction(4); D26=Fraction(26); SO=Fraction(10)
v = SO*D4*D4 - D26 + K + P + F
val = float(v); target = 137.036
print(f"S475 COMPLETE. Fine structure 1/alpha = {val:.4f}; target {target}; match {abs(val-target)/target*100:.4f}%.")
