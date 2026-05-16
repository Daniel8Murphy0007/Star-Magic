from fractions import Fraction
D26=Fraction(26); D4=Fraction(4); F=Fraction(1,10); P=Fraction(5,6); K=Fraction(25,12)
# Saturn orbital period 29.46 yr
v = D26 + D4 - F*K - 3*F*P - F*F*K - F*F*P - F*F - F*F*F
val = float(v); target = 29.46
print(f"S522 COMPLETE. Saturn period = {val:.6f}; closure = D_crit+D_phys-F*K-3*F*Phi_res-F^2*K-F^2*Phi_res-F^2-F^3; target {target}; match {abs(val-target)/target*100:.4f}%.")
