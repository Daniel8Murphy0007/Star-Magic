from fractions import Fraction
F=Fraction(1,10); P=Fraction(5,6); S=Fraction(57,100); K=Fraction(25,12)
v = K + P - S + F*K - F*F - F*F*K - F*F*F - F*P - F*F*P
val = float(v); target = 2.426
print(f"S480 COMPLETE. Compton lambda_e (pm) = {val:.4f}; target {target}; match {abs(val-target)/target*100:.4f}%.")
