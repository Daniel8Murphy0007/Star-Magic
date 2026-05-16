from fractions import Fraction
F=Fraction(1,10); Ph=Fraction(5,6); K=Fraction(25,12)
v = Ph + F*K + F*Ph + F*F*K + F*F*Ph; target = Fraction(1161,1000)
print(f"S523 COMPLETE. Pareto 80/20 exponent = {float(v):.6f}; closure = Phi_res+F*K+F*Phi_res+F^2*K+F^2*Phi_res; target {float(target)}; match {abs(float(v)-float(target))/float(target)*100:.4f}%.")
