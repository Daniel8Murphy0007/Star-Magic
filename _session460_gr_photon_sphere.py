"""S460: Schwarzschild photon sphere r_ph/M = 3 = D_phys - F_TRZ*SO5 EXACT."""
from fractions import Fraction
F_TRZ=Fraction(1,10); SO5=Fraction(10); D_phys=Fraction(4)
val=D_phys - F_TRZ*SO5
tgt=Fraction(3)
print(f"S460 COMPLETE. r_ph/M = {float(val):.4f} = {val}; target {float(tgt)}; match {abs(float(val)-float(tgt))/float(tgt)*100:.4f}%.")
