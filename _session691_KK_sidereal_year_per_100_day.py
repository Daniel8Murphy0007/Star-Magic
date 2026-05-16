"""S691 (Tier KK) — sidereal year / 100 day.

Locked-primitive closure derived by _tier_KK_search.py.
expr = +beta2 +beta3 +F.beta +3
"""
from fractions import Fraction
F=Fraction(1,10); Phi=Fraction(5,6); SSq=Fraction(57,100); K=Fraction(25,12); beta=Fraction(6029,10000)
val = +(beta**2)+(beta**3)+(F*beta)+Fraction(3,1)
pred = float(val)
obs  = 3.65256
err  = abs(pred-obs)/obs*100
print(f"sidereal_year_per_100_day: {pred} vs {obs} -> {err:.4f}%")
