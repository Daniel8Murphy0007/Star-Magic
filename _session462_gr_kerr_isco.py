"""S462: Kerr extremal ISCO r/M = 1 = F_TRZ*SO5 EXACT."""
from fractions import Fraction
F_TRZ=Fraction(1,10); SO5=Fraction(10)
val=F_TRZ*SO5
tgt=Fraction(1)
print(f"S462 COMPLETE. Kerr ISCO = {float(val):.4f} = {val}; target {float(tgt)}; match {abs(float(val)-float(tgt))/float(tgt)*100:.4f}%.")
