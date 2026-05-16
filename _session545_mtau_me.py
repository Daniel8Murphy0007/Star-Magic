from fractions import Fraction
F=Fraction(1,10); Ph=Fraction(5,6); K=Fraction(25,12); SSq=Fraction(57,100)
Dphys=4; Dbsfg=6; Dcrit=26; Nch=9; SO5=10; A5=60
v = A5*A5 - A5*Dphys + A5*F*Dphys + A5*Ph + Nch*Dphys + F*SO5*Dphys + Nch*F*Dphys
target = 3477.23
print(f"S545 COMPLETE. m_tau/m_e = {float(v):.4f}; closure = A_5^2-A_5*D_phys+A_5*F*D_phys+A_5*Phi_res+N_ch*D_phys+F*SO5*D_phys+N_ch*F*D_phys; target {target:.4f}; match {abs(float(v)-target)/target*100:.4f}%.")
