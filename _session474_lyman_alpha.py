from fractions import Fraction
F=Fraction(1,10); P=Fraction(5,6); K=Fraction(25,12); SO=Fraction(10)
v = SO*F + F*K - F*F*F + F*F*K - F*F*P
val = float(v); target = 1.216
print(f"S474 COMPLETE. Lyman-alpha (x100 nm) = {val:.4f}; target {target}; match {abs(val-target)/target*100:.4f}%.")
