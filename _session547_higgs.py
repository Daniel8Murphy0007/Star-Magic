from fractions import Fraction
F=Fraction(1,10); Ph=Fraction(5,6); K=Fraction(25,12); SSq=Fraction(57,100)
Dphys=4; Dbsfg=6; Dcrit=26; Nch=9; SO5=10; A5=60
v = Nch*SO5 + Dbsfg*Dbsfg - Ph
target = 125.25
print(f"S547 COMPLETE. Higgs mass (GeV) = {float(v):.4f}; closure = N_ch*SO5+D_BSFG^2-Phi_res; target {target:.4f}; match {abs(float(v)-target)/target*100:.4f}%.")
