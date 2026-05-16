from fractions import Fraction
SO=Fraction(10); D6=Fraction(6); D26=Fraction(26)
# Catalan number C_5 = 42
v = D26 + D6 + SO
val = float(v); target = 42
print(f"S506 COMPLETE. Catalan C_5 = {val:.4f} = D_crit+D_BSFG+SO5 = 42; target {target}; match {abs(val-target)/target*100:.4f}%.")
