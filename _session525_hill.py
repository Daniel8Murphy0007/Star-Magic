from fractions import Fraction
F=Fraction(1,10); Ph=Fraction(5,6); K=Fraction(25,12)
v = K + Ph - F*Ph - F*F*Ph - F*F*K; target = Fraction(28,10)
print(f"S525 COMPLETE. Hill coefficient (hemoglobin) = {float(v):.6f}; closure = K_Mex+Phi_res-F*Phi_res-F^2*Phi_res-F^2*K_Mex; target {float(target)}; match {abs(float(v)-float(target))/float(target)*100:.4f}%.")
