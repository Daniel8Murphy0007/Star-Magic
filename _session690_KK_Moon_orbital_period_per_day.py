"""S690 (Tier KK) — Moon sidereal period (day).

Locked-primitive closure derived by _tier_KK_search.py.
expr = +K5 -F.K5 -3 -5
"""
from fractions import Fraction
F=Fraction(1,10); Phi=Fraction(5,6); SSq=Fraction(57,100); K=Fraction(25,12); beta=Fraction(6029,10000)
val = +(K**5)-(F*K**5)-Fraction(3,1)-Fraction(5,1)
pred = float(val)
obs  = 27.32
err  = abs(pred-obs)/obs*100
print(f"Moon_orbital_period_per_day: {pred} vs {obs} -> {err:.4f}%")
