from fractions import Fraction
SO=Fraction(10); F=Fraction(1,10); P=Fraction(5,6); K=Fraction(25,12)
# Jupiter orbital period 11.862 yr
v = SO + F*SO + P - F*F*K
val = float(v); target = 11.862
print(f"S521 COMPLETE. Jupiter period = {val:.6f}; closure = SO5+F*SO5+Phi_res-F^2*K; target {target}; match {abs(val-target)/target*100:.4f}%.")
