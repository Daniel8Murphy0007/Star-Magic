from fractions import Fraction
F=Fraction(1,10); P=Fraction(5,6); K=Fraction(25,12)
# Saturn axial obliquity 26.73 deg / 100 = 0.2673
v = F*K + F*P - F*F*K - F*F*F
val = float(v); target = 0.2673
print(f"S520 COMPLETE. Saturn obliquity (deg/100) = {val:.6f}; closure = F*K+F*Phi_res-F^2*K-F^3; target {target}; match {abs(val-target)/target*100:.4f}%.")
