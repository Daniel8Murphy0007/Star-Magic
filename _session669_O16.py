from fractions import Fraction
F=Fraction(1,10); K=Fraction(25,12); beta=Fraction(6029,10000)
v = F*K**4 + F*K**5 + beta**4 + F*beta**2 + 2
print(f"O16 BE/A: {float(v):.6f} vs 7.976200 -> {abs(float(v)-7.9762)/7.9762*100:.4f}%")
