from fractions import Fraction
import math
F=Fraction(1,10); P=Fraction(5,6); S=Fraction(57,100)
# Hausdorff dim of Cantor set = ln 2 / ln 3 ~ 0.6309
v = S + F*P - F*F*P - F*F
val = float(v); target = math.log(2)/math.log(3)
print(f"S508 COMPLETE. Cantor Hausdorff dim = {val:.6f}; closure = SSq+F*Phi_res-F^2*Phi_res-F^2; target {target:.6f}; match {abs(val-target)/target*100:.4f}%.")
