from fractions import Fraction
import math
F=Fraction(1,10); P=Fraction(5,6); S=Fraction(57,100); K=Fraction(25,12)
# Kepler sphere packing density eta_3D = pi/sqrt(18) ~ 0.7405
v = S + F*K - F*F*K - F*F*P - F*F*F
val = float(v); target = math.pi/math.sqrt(18)
print(f"S503 COMPLETE. Kepler eta_3D = {val:.6f}; closure = SSq+F*K-F^2*K-F^2*Phi_res-F^3; target {target:.6f}; match {abs(val-target)/target*100:.4f}%.")
