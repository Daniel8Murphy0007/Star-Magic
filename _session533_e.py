from fractions import Fraction
import math
F=Fraction(1,10); Ph=Fraction(5,6); K=Fraction(25,12)
v = K + Ph - F*K + F*F*K - F*F*Ph
target = math.e
print(f"S533 COMPLETE. e = {float(v):.6f}; closure = K_Mex+Phi_res-F*K+F^2*K-F^2*Phi_res; target {target:.6f}; match {abs(float(v)-target)/target*100:.4f}%.")
