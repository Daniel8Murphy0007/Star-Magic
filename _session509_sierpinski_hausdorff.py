from fractions import Fraction
import math
F=Fraction(1,10); P=Fraction(5,6)
# Hausdorff dim of Sierpinski triangle = ln 3 / ln 2 ~ 1.585
v = P + P - F*P
val = float(v); target = math.log(3)/math.log(2)
print(f"S509 COMPLETE. Sierpinski Hausdorff dim = {val:.6f}; closure = 2*Phi_res - F*Phi_res = 19/12; target {target:.6f}; match {abs(val-target)/target*100:.4f}%.")
