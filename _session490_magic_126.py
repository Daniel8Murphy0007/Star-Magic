from fractions import Fraction
D26=Fraction(26); SO=Fraction(10)
v = D26 + SO*SO
val = float(v); target = 126
print(f"S490 COMPLETE. Magic number n-drip = {val:.4f} = 126 = D_crit + SO5^2; target {target}; match {abs(val-target)/target*100:.4f}%.")
