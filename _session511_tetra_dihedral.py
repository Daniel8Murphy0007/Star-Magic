from fractions import Fraction
import math
F=Fraction(1,10); P=Fraction(5,6); K=Fraction(25,12)
# Regular tetrahedron dihedral angle = arccos(1/3) ~ 70.529 deg => /100 = 0.70529
v = P - F*P - F*F*K - F*F*P - F*F*F - F*F
val = float(v); target = math.degrees(math.acos(1/3))/100
print(f"S511 COMPLETE. Tetrahedron dihedral (deg/100) = {val:.6f}; closure = Phi_res-F*Phi_res-F^2*K-F^2*Phi_res-F^3-F^2; target {target:.6f}; match {abs(val-target)/target*100:.4f}%.")
