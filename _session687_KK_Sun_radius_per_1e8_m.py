"""S687 (Tier KK) — R_sun / 10^8 m.

Locked-primitive closure derived by _tier_KK_search.py.
expr = -F2.beta2 -F2.beta3 +2 +5
"""
from fractions import Fraction
F=Fraction(1,10); Phi=Fraction(5,6); SSq=Fraction(57,100); K=Fraction(25,12); beta=Fraction(6029,10000)
val = -(F**2*beta**2)-(F**2*beta**3)+Fraction(2,1)+Fraction(5,1)
pred = float(val)
obs  = 6.96
err  = abs(pred-obs)/obs*100
print(f"Sun_radius_per_1e8_m: {pred} vs {obs} -> {err:.4f}%")
