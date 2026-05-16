"""S396 von Klitzing constant R_K = h/e^2 = 25812.807 Ohm; log10 = 4.4118"""
import math
F_TRZ=0.1; Phi_res=5/6; SSq=0.57; D_phys=4
val = D_phys + SSq*Phi_res - SSq*F_TRZ + F_TRZ*F_TRZ*Phi_res
obs = math.log10(25812.807)
print(f"S396 COMPLETE. log10(R_K) = D_phys+SSq*Phi_res-SSq*F_TRZ+F_TRZ^2*Phi_res = {val:.4f}; obs (log10 h/e^2) = {obs:.4f}; match {abs(val-obs)/obs*100:.3f}%.")
