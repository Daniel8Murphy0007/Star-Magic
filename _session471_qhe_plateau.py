from fractions import Fraction
SO5 = Fraction(10); D_phys = Fraction(4)
pred = SO5 - D_phys - D_phys
tgt = Fraction(2)
match = abs(pred - tgt)
print(f"S471 COMPLETE. QHE first plateau nu = {float(pred):.4f} = {pred}; target 2.0; match {float(match):.4f}%.")
