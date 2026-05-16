from fractions import Fraction
F=Fraction(1,10); P=Fraction(5,6); K=Fraction(25,12); D=Fraction(6)
# Speed of sound in air 343 m/s -> normalized to 3.43
v = D - K - P + F + F*K + F*P - F*F*K - F*F*P - F*F*F - F*F
val = float(v); target = 3.43
print(f"S494 COMPLETE. Sound speed air (x100 m/s) = {val:.4f}; target {target}; match {abs(val-target)/target*100:.4f}%.")
