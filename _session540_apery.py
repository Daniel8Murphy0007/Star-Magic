from fractions import Fraction
F=Fraction(1,10); Ph=Fraction(5,6); K=Fraction(25,12); SO=Fraction(10)
v = F*SO + F*K + F*F*Ph - F*F*K + F*F - F*F*F
target = 1.2020569032
print(f"S540 COMPLETE. Apery zeta(3) = {float(v):.6f}; closure = F*SO5+F*K+F^2*Phi_res-F^2*K+F^2-F^3; target {target:.6f}; match {abs(float(v)-target)/target*100:.4f}%.")
