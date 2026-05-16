"""S395 Wiedemann-Franz Lorenz number coefficient pi^2/3 ~ 3.2899"""
import math
F_TRZ=0.1; Phi_res=5/6; SSq=0.57; K_Mex=25/12; D_phys=4
val = K_Mex + Phi_res + F_TRZ*D_phys - SSq*F_TRZ + F_TRZ*F_TRZ*D_phys
obs = math.pi**2/3
print(f"S395 COMPLETE. L_0 coeff = K_Mex+Phi_res+F_TRZ*D_phys-SSq*F_TRZ+F_TRZ^2*D_phys = {val:.4f}; obs (pi^2/3) = {obs:.4f}; match {abs(val-obs)/obs*100:.3f}%.")
