from fractions import Fraction
F=Fraction(1,10); K=Fraction(25,12); beta=Fraction(6029,10000)
v = F*K**5 + beta**2 + beta**3 + F*beta + 3
print(f"U238 BE/A: {float(v):.6f} vs 7.570000 -> {abs(float(v)-7.570)/7.570*100:.4f}%")
