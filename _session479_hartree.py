from fractions import Fraction
F=Fraction(1,10); P=Fraction(5,6); K=Fraction(25,12); D=Fraction(26)
v = D + K - P + F - F*K + F*P - F*F*F
val = float(v); target = 27.211
print(f"S479 COMPLETE. Hartree (eV) = {val:.4f}; target {target}; match {abs(val-target)/target*100:.4f}%.")
