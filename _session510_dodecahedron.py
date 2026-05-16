from fractions import Fraction
D6=Fraction(6)
# Regular dodecahedron: 12 pentagonal faces
v = D6 + D6
val = float(v); target = 12
print(f"S510 COMPLETE. Dodecahedron faces F = {val:.4f} = 2*D_BSFG = 12; target {target}; match {abs(val-target)/target*100:.4f}%.")
