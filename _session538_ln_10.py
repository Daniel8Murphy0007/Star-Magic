from fractions import Fraction
import math
F=Fraction(1,10); K=Fraction(25,12)
v = K + F*K + F*F + F*F*F
target = math.log(10)
print(f"S538 COMPLETE. ln 10 = {float(v):.6f}; closure = K_Mex+F*K_Mex+F^2+F^3; target {target:.6f}; match {abs(float(v)-target)/target*100:.4f}%.")
