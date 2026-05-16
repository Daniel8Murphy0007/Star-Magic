"""S578: Earth-Sun distance 149.6 Gm = D_c*D_B - D_phys - K_Mex - F_TRZ*D_phys + F_TRZ*Phi_res (EXACT)
       156 - 4 - 25/12 - 0.4 + 1/12 = 152 - 24/12 - 0.4 = 152 - 2 - 0.4 = 149.6"""
from fractions import Fraction
F=Fraction(1,10); Ph=Fraction(5,6); K=Fraction(25,12); D_p=4; D_B=6; D_c=26
v=D_c*D_B-D_p-K-F*D_p+F*Ph
target=Fraction(1496,10)
print(f"AU (Gm): {float(v):.4f} vs {float(target)} -> {abs(float(v-target))/float(target)*100:.4f}% [EXACT]" if v==target else f"{float(v)}")
