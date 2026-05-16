"""S581: Sidereal year days = 365.25 = N_ch*A_5 - D_phys*A_5 + A_5 + D_phys + K_Mex - Phi_res (EXACT)
       = 540 - 240 + 60 + 4 + 25/12 - 10/12 = 364 + 15/12 = 364 + 1.25 = 365.25"""
from fractions import Fraction
Ph=Fraction(5,6); K=Fraction(25,12); D_p=4; N=9; A=60
v=N*A-D_p*A+A+D_p+K-Ph
target=Fraction(36525,100)
print(f"Year days: {float(v)} vs {float(target)} -> {abs(float(v-target))/float(target)*100:.4f}% [EXACT]" if v==target else f"{float(v)}")
