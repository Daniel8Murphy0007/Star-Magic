"""S689 (Tier KK) — R_earth / 10^6 m.

Locked-primitive closure derived by _tier_KK_search.py.
expr = +beta2 -2 +3 +5
"""
from fractions import Fraction
F=Fraction(1,10); Phi=Fraction(5,6); SSq=Fraction(57,100); K=Fraction(25,12); beta=Fraction(6029,10000)
val = +(beta**2)-Fraction(2,1)+Fraction(3,1)+Fraction(5,1)
pred = float(val)
obs  = 6.371
err  = abs(pred-obs)/obs*100
print(f"Earth_radius_per_1e6_m: {pred} vs {obs} -> {err:.4f}%")
