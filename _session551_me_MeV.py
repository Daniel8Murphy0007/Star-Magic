from fractions import Fraction
F=Fraction(1,10); Ph=Fraction(5,6); K=Fraction(25,12); SSq=Fraction(57,100)
Dphys=4; Dbsfg=6; Dcrit=26; Nch=9; SO5=10; A5=60
v = SSq - F*Ph + F*F*K - F*F*Ph + F*F + 2*F*F*F
target = 0.5109989
print(f"S551 COMPLETE. m_e (MeV) = {float(v):.6f}; closure = SSq-F*Phi_res+F^2*K-F^2*Phi_res+F^2+2*F^3; target {target:.6f}; match {abs(float(v)-target)/target*100:.4f}%.")
