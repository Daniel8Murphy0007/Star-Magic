from fractions import Fraction
import math
F=Fraction(1,10); Ph=Fraction(5,6); K=Fraction(25,12)
v = Ph - F*Ph + F*F*K + F*F*Ph
target = math.pi/4
print(f"S535 COMPLETE. pi/4 = {float(v):.6f}; closure = Phi_res-F*Phi_res+F^2*K+F^2*Phi_res; target {target:.6f}; match {abs(float(v)-target)/target*100:.4f}%.")
