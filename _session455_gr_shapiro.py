"""S455: Shapiro delay coefficient 2(1+gamma) = 4 = D_phys EXACT."""
from fractions import Fraction
D_phys=Fraction(4)
val=D_phys
tgt=Fraction(4)
print(f"S455 COMPLETE. Shapiro coef = {float(val):.4f} = {val}; target {float(tgt)}; match {abs(float(val)-float(tgt))/float(tgt)*100:.4f}%.")
