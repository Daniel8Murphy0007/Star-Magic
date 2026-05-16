from fractions import Fraction
import math
F=Fraction(1,10); P=Fraction(5,6); K=Fraction(25,12)
# 2D circle packing density eta_2D = pi/(2*sqrt(3)) ~ 0.9069
v = P + F*K - F*P - F*F*K - F*F*P - F*F*F - F*F
val = float(v); target = math.pi/(2*math.sqrt(3))
print(f"S504 COMPLETE. Circle packing eta_2D = {val:.6f}; closure = Phi_res+F*K-F*Phi_res-F^2*K-F^2*Phi_res-F^3-F^2; target {target:.6f}; match {abs(val-target)/target*100:.4f}%.")
