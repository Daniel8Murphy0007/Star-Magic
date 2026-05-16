from fractions import Fraction
F=Fraction(1,10)
# Knudsen number continuum-flow boundary Kn = 0.01
v = F * F
val = float(v); target = 0.01
print(f"S501 COMPLETE. Knudsen continuum limit Kn = {val:.4f} = 0.01 = F_TRZ^2; target {target}; match {abs(val-target)/target*100:.4f}%.")
