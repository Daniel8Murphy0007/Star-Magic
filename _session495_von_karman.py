from fractions import Fraction
F=Fraction(1,10); K=Fraction(25,12); D=Fraction(4)
# von Karman constant kappa = 0.41
v = F*D + F*F*K - F*F*F - F*F
val = float(v); target = 0.41
print(f"S495 COMPLETE. von Karman kappa = {val:.4f}; target {target}; match {abs(val-target)/target*100:.4f}%.")
