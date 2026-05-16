from fractions import Fraction
F=Fraction(1,10); Ph=Fraction(5,6); K=Fraction(25,12); SSq=Fraction(57,100)
Dphys=4; Dbsfg=6; Dcrit=26; Nch=9; SO5=10; A5=60
v = F*SO5 + F*F*F + F*F*F*F
target = 1.001378
print(f"S549 COMPLETE. m_n/m_p = {float(v):.6f}; closure = F*SO5+F^3+F^4; target {target:.6f}; match {abs(float(v)-target)/target*100:.4f}%.")
