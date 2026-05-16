from fractions import Fraction
A5=Fraction(60); D26=Fraction(26); D4=Fraction(4)
v = A5 + D26 - D4
val = float(v); target = 82
print(f"S489 COMPLETE. Magic number Pb = {val:.4f} = 82 = A_5 + D_crit - D_phys; target {target}; match {abs(val-target)/target*100:.4f}%.")
