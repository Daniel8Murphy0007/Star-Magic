from fractions import Fraction
F=Fraction(1,10); Ph=Fraction(5,6); K=Fraction(25,12); SSq=Fraction(57,100)
Dphys=4; Dbsfg=6; Dcrit=26; Nch=9; SO5=10; A5=60
v = Nch*SO5*SO5 + Nch*Dphys + K + 2*F*Ph
target = 938.27208816
print(f"S550 COMPLETE. m_p (MeV) = {float(v):.4f}; closure = N_ch*SO5^2+N_ch*D_phys+K_Mex+2*F*Phi_res; target {target:.4f}; match {abs(float(v)-target)/target*100:.4f}%.")
