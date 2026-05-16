from fractions import Fraction
F=Fraction(1,10); P=Fraction(5,6); K=Fraction(25,12); S=Fraction(57,100)
# Mars Bode-law distance ~1.6 AU
v = K - S + F*K - F*P - F*F*K - F*F*P
val = float(v); target = 1.6
print(f"S518 COMPLETE. Mars Bode = {val:.6f}; closure = K_Mex-SSq+F*K-F*Phi_res-F^2*K-F^2*Phi_res; target {target}; match {abs(val-target)/target*100:.4f}%.")
