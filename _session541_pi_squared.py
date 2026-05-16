from fractions import Fraction
import math
F=Fraction(1,10); Ph=Fraction(5,6); K=Fraction(25,12); SO=Fraction(10)
v = SO - F - F*F*K - F*F*Ph
target = math.pi**2
print(f"S541 COMPLETE. pi^2 = {float(v):.6f}; closure = SO5-F_TRZ-F^2*K-F^2*Phi_res; target {target:.6f}; match {abs(float(v)-target)/target*100:.4f}%.")
