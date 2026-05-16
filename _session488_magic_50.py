from fractions import Fraction
A5=Fraction(60); SO=Fraction(10)
v = A5 - SO
val = float(v); target = 50
print(f"S488 COMPLETE. Magic number Sn = {val:.4f} = 50 = A_5 - SO5; target {target}; match {abs(val-target)/target*100:.4f}%.")
