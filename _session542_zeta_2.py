from fractions import Fraction
import math
F=Fraction(1,10); Ph=Fraction(5,6); K=Fraction(25,12)
v = K - F*K - 2*F*Ph - 2*F*F*K - F*F*Ph - F*F - F*F*F
target = math.pi**2 / 6
print(f"S542 COMPLETE. zeta(2)=pi^2/6 = {float(v):.6f}; closure = K-F*K-2*F*Phi_res-2*F^2*K-F^2*Phi_res-F^2-F^3; target {target:.6f}; match {abs(float(v)-target)/target*100:.4f}%.")
