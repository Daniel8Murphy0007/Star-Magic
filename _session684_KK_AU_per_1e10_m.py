"""S684 (Tier KK) — 1 AU / 10^10 m.

Locked-primitive closure derived by _tier_KK_search.py.
expr = +K4 -beta -beta3 -3
"""
from fractions import Fraction
F=Fraction(1,10); Phi=Fraction(5,6); SSq=Fraction(57,100); K=Fraction(25,12); beta=Fraction(6029,10000)
val = +(K**4)-(beta)-(beta**3)-Fraction(3,1)
pred = float(val)
obs  = 14.96
err  = abs(pred-obs)/obs*100
print(f"AU_per_1e10_m: {pred} vs {obs} -> {err:.4f}%")
