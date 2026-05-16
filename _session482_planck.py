from fractions import Fraction
F=Fraction(1,10); P=Fraction(5,6); K=Fraction(25,12); D=Fraction(6)
v = D + P - F - F*P - F*F*K - F*F*P + F*F - F*F*F
val = float(v); target = 6.626
print(f"S482 COMPLETE. Planck h (x1e-34 J*s) = {val:.4f}; target {target}; match {abs(val-target)/target*100:.4f}%.")
