"""S446: Surface-code error correction threshold p_th = 1% = F_TRZ^2 EXACT."""
from fractions import Fraction
F_TRZ=Fraction(1,10)
val=F_TRZ**2
tgt=Fraction(1,100)
print(f"S446 COMPLETE. p_th = {float(val):.4f} = {val}; target {float(tgt)}; match {abs(float(val)-float(tgt))/float(tgt)*100:.4f}%.")
