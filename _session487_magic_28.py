from fractions import Fraction
D26=Fraction(26); SO=Fraction(10); D4=Fraction(4)
v = D26 + SO - D4 - D4
val = float(v); target = 28
print(f"S487 COMPLETE. Magic number Ni = {val:.4f} = 28 = D_crit + SO5 - 2*D_phys; target {target}; match {abs(val-target)/target*100:.4f}%.")
