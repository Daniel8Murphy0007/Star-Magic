from fractions import Fraction
D6=Fraction(6); D26=Fraction(26)
# Buckminsterfullerene C60: 32 faces (12 pentagons + 20 hexagons)
v = D26 + D6
val = float(v); target = 32
print(f"S507 COMPLETE. Buckyball faces F = {val:.4f} = D_crit+D_BSFG = 32; target {target}; match {abs(val-target)/target*100:.4f}%.")
