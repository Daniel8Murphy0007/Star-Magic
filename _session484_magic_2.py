from fractions import Fraction
SO=Fraction(10); D=Fraction(4)
v = SO - D - D
val = float(v); target = 2
print(f"S484 COMPLETE. Magic number He = {val:.4f} = 2; target {target}; match {abs(val-target)/target*100:.4f}%.")
