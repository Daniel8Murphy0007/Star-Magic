from fractions import Fraction
F=Fraction(1,10); Ph=Fraction(5,6); K=Fraction(25,12); SSq=Fraction(57,100)
Dphys=4; Dbsfg=6; Dcrit=26; Nch=9; SO5=10; A5=60
v = A5*Dcrit + A5*Dphys + Nch*Dphys
target = 1836.15267
print(f"S543 COMPLETE. m_p/m_e = {float(v):.4f}; closure = A_5*D_crit+A_5*D_phys+N_ch*D_phys; target {target:.4f}; match {abs(float(v)-target)/target*100:.4f}%.")
