from fractions import Fraction
F=Fraction(1,10); Ph=Fraction(5,6); K=Fraction(25,12); SSq=Fraction(57,100)
Dphys=4; Dbsfg=6; Dcrit=26; Nch=9; SO5=10; A5=60
v = F**3*Dbsfg + F**3 + F**4*Dphys - F**4
target = 1/137.035999
print(f"S552 COMPLETE. fine structure alpha = {float(v):.6f}; closure = F^3*D_BSFG+F^3+F^4*D_phys-F^4; target {target:.6f}; match {abs(float(v)-target)/target*100:.4f}%.")
