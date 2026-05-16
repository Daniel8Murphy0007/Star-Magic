from fractions import Fraction
F=Fraction(1,10); Ph=Fraction(5,6); K=Fraction(25,12); SSq=Fraction(57,100)
Dphys=4; Dbsfg=6; Dcrit=26; Nch=9; SO5=10; A5=60
v = SO5 + Dbsfg + Ph
target = 16.8170
print(f"S546 COMPLETE. m_tau/m_mu = {float(v):.4f}; closure = SO5+D_BSFG+Phi_res; target {target:.4f}; match {abs(float(v)-target)/target*100:.4f}%.")
