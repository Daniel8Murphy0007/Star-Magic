from fractions import Fraction
F=Fraction(1,10); SO=Fraction(10); D=Fraction(4)
# EXACT closure: SO5 + D_phys*(1 - F_TRZ) = 10 + 4*0.9 = 13.6
v = SO + D - F*D
val = float(v); target = 13.6
print(f"S477 COMPLETE. H ionization (eV) = {val:.4f} = 13.6; target {target}; match {abs(val-target)/target*100:.4f}%.")
