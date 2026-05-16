"""S685 (Tier KK) — M_sun / 10^29 kg.

Locked-primitive closure derived by _tier_KK_search.py.
expr = +K4 -2 +3
"""
from fractions import Fraction
F=Fraction(1,10); Phi=Fraction(5,6); SSq=Fraction(57,100); K=Fraction(25,12); beta=Fraction(6029,10000)
val = +(K**4)-Fraction(2,1)+Fraction(3,1)
pred = float(val)
obs  = 19.89
err  = abs(pred-obs)/obs*100
print(f"Sun_mass_per_1e29_kg: {pred} vs {obs} -> {err:.4f}%")
