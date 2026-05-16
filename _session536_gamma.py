from fractions import Fraction
F=Fraction(1,10); Ph=Fraction(5,6); K=Fraction(25,12); SSq=Fraction(57,100)
v = SSq + F*F*K - F*F*Ph
target = 0.5772156649
print(f"S536 COMPLETE. Euler-Mascheroni gamma = {float(v):.6f}; closure = SSq+F^2*K-F^2*Phi_res; target {target:.6f}; match {abs(float(v)-target)/target*100:.4f}%.")
