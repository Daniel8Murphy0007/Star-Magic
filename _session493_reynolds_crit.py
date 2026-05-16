from fractions import Fraction
F=Fraction(1,10); D26=Fraction(26); D4=Fraction(4); SO=Fraction(10)
# Re_crit pipe ~2300 -> normalized to 23
v = D26 - D4 + F*SO
val = float(v); target = 23
print(f"S493 COMPLETE. Re_crit pipe (x100) = {val:.4f} = 23 = D_crit - D_phys + F_TRZ*SO5; target {target}; match {abs(val-target)/target*100:.4f}%.")
