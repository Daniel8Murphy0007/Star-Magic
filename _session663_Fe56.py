from fractions import Fraction
F=Fraction(1,10); K=Fraction(25,12); beta=Fraction(6029,10000)
v = F*K**5 - beta**4 + 2 + 3
print(f"Fe56 BE/A: {float(v):.6f} vs 8.790300 -> {abs(float(v)-8.7903)/8.7903*100:.4f}%")
