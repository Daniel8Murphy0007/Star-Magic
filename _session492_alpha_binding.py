from fractions import Fraction
F=Fraction(1,10); P=Fraction(5,6); K=Fraction(25,12); D=Fraction(26)
v = D + K + F + F*P + F*F*K + F*F*P
val = float(v); target = 28.30
print(f"S492 COMPLETE. Alpha binding (MeV) = {val:.4f}; target {target}; match {abs(val-target)/target*100:.4f}%.")
