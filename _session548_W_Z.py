from fractions import Fraction
F=Fraction(1,10); Ph=Fraction(5,6); K=Fraction(25,12); SSq=Fraction(57,100)
Dphys=4; Dbsfg=6; Dcrit=26; Nch=9; SO5=10; A5=60
v = SSq + F*K + F*Ph + F*F*K
target = 0.88153
print(f"S548 COMPLETE. M_W/M_Z = {float(v):.6f}; closure = SSq+F*K+F*Phi_res+F^2*K; target {target:.6f}; match {abs(float(v)-target)/target*100:.4f}%.")
