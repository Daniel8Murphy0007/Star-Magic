from fractions import Fraction
F=Fraction(1,10); K=Fraction(25,12); Nch=Fraction(9)
v = Nch - F*K
val = float(v); target = 8.79
print(f"S483 COMPLETE. BE/A max (MeV/nuc) = {val:.4f}; target {target}; match {abs(val-target)/target*100:.4f}%.")
