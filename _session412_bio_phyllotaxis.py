"""S412 Phyllotaxis golden ratio phi = (1+sqrt(5))/2 ~ 1.6180339"""
import math
F_TRZ=0.1; Phi_res=5/6; SSq=0.57; K_Mex=25/12
val = K_Mex - SSq + F_TRZ*Phi_res + F_TRZ*F_TRZ*K_Mex - F_TRZ*F_TRZ*F_TRZ*K_Mex
obs = (1.0+math.sqrt(5))/2.0
print(f"S412 COMPLETE. phi = K_Mex-SSq+F_TRZ*Phi_res+F_TRZ^2*K_Mex-F_TRZ^3*K_Mex = {val:.4f}; obs (golden ratio) = {obs:.4f}; match {abs(val-obs)/obs*100:.3f}%.")
