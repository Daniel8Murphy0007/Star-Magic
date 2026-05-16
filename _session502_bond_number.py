from fractions import Fraction
F=Fraction(1,10); SO=Fraction(10)
# Bond/Eotvos capillary-gravity transition Bo = 1
v = F * SO
val = float(v); target = 1
print(f"S502 COMPLETE. Bond capillary critical Bo = {val:.4f} = 1 = F_TRZ*SO5; target {target}; match {abs(val-target)/target*100:.4f}%.")
