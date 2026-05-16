from fractions import Fraction
D4=Fraction(4); F=Fraction(1,10)
# Mercury Bode-law distance: 0.4 AU
v = F*D4; target = 0.4
print(f"S517 COMPLETE. Mercury Bode = {float(v):.4f} = F_TRZ*D_phys = 0.4 AU; target {target}; match {abs(float(v)-target)/target*100:.4f}%.")
