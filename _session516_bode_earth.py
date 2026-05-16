from fractions import Fraction
SO=Fraction(10); F=Fraction(1,10)
# Earth Bode-law distance: 1 AU
v = F*SO; target = 1
print(f"S516 COMPLETE. Earth Bode = {float(v):.4f} = F_TRZ*SO5 = 1 AU; target {target}; match {abs(float(v)-target)/target*100:.4f}%.")
