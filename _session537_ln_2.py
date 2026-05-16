from fractions import Fraction
import math
F=Fraction(1,10); Ph=Fraction(5,6); K=Fraction(25,12)
v = 2*F + Ph - F*K - F*Ph - F*F*K - 2*F*F*Ph - F*F*F - F*F
target = math.log(2)
print(f"S537 COMPLETE. ln 2 = {float(v):.6f}; closure = 2F+Phi_res-F*K-F*Phi_res-F^2*K-2*F^2*Phi_res-F^3-F^2; target {target:.6f}; match {abs(float(v)-target)/target*100:.4f}%.")
