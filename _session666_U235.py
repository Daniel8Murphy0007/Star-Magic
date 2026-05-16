from fractions import Fraction
F=Fraction(1,10); K=Fraction(25,12); beta=Fraction(6029,10000)
v = F*K**5 + beta + F*beta + 3
print(f"U235 BE/A: {float(v):.6f} vs 7.591000 -> {abs(float(v)-7.591)/7.591*100:.4f}%")
