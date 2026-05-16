from fractions import Fraction
import math
F=Fraction(1,10); Ph=Fraction(5,6); K=Fraction(25,12); SO=Fraction(10); Db=Fraction(6)
v = Db + K - F*SO + F*Ph + F*K + F*F*K
target = math.e**2
print(f"S534 COMPLETE. e^2 = {float(v):.6f}; closure = D_BSFG+K-F*SO5+F*Phi_res+F*K+F^2*K; target {target:.6f}; match {abs(float(v)-target)/target*100:.4f}%.")
