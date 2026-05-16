"""S686 (Tier KK) — Earth orbital v (km/s).

Locked-primitive closure derived by _tier_KK_search.py.
expr = +K5 -F.K5 -beta -5
"""
from fractions import Fraction
F=Fraction(1,10); Phi=Fraction(5,6); SSq=Fraction(57,100); K=Fraction(25,12); beta=Fraction(6029,10000)
val = +(K**5)-(F*K**5)-(beta)-Fraction(5,1)
pred = float(val)
obs  = 29.78
err  = abs(pred-obs)/obs*100
print(f"Earth_orbit_v_per_km_s: {pred} vs {obs} -> {err:.4f}%")
