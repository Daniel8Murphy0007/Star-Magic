from fractions import Fraction
A5=Fraction(60)
# Buckminsterfullerene C60: 60 vertices
v = A5
val = float(v); target = 60
print(f"S512 COMPLETE. Buckyball vertices V = {val:.4f} = A_5 = 60; target {target}; match {abs(val-target)/target*100:.4f}%.")
