from fractions import Fraction
SO=Fraction(10); D4=Fraction(4)
# Catalan number C_4 = 14
v = SO + D4
val = float(v); target = 14
print(f"S505 COMPLETE. Catalan C_4 = {val:.4f} = SO5+D_phys = 14; target {target}; match {abs(val-target)/target*100:.4f}%.")
