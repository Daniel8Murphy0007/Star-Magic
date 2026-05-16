from fractions import Fraction
F=Fraction(1,10); K=Fraction(25,12); beta=Fraction(6029,10000)
v = F*K**5 + beta + beta**4 + F*beta**3 + 3
print(f"C12 BE/A: {float(v):.6f} vs 7.680200 -> {abs(float(v)-7.6802)/7.6802*100:.4f}%")
