from fractions import Fraction
F=Fraction(1,10); Ph=Fraction(5,6); K=Fraction(25,12); SSq=Fraction(57,100)
v = SSq - F*K + F*F*K + F*F*Ph - F*F*F - F*F
# Phyllotaxis golden angle 137.5/360 = 0.38194 ≈ 1/phi^2
target = 137.5/360
print(f"S527 COMPLETE. Golden phyllotaxis (137.5/360) = {float(v):.6f}; closure = SSq-F*K+F^2*K+F^2*Phi_res-F^3-F^2; target {target:.6f}; match {abs(float(v)-target)/target*100:.4f}%.")
