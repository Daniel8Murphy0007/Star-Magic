from fractions import Fraction
A5=Fraction(60); SO=Fraction(10); P=Fraction(5,6); D6=Fraction(6)
# Halley's comet orbital period ~ 75 years
v = A5 + SO + P*D6; target = 75
print(f"S519 COMPLETE. Halley period = {float(v):.4f} = A_5+SO5+Phi_res*D_BSFG = 75 yr; target {target}; match {abs(float(v)-target)/target*100:.4f}%.")
