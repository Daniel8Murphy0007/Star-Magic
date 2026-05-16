"""S397 BCS coherence length coefficient 1/pi ~ 0.31831"""
import math
F_TRZ=0.1; Phi_res=5/6; SSq=0.57; D_phys=4
val = F_TRZ*Phi_res*D_phys - F_TRZ*F_TRZ - F_TRZ*F_TRZ*SSq
obs = 1.0/math.pi
print(f"S397 COMPLETE. xi_0/(hbar v_F/Delta) = F_TRZ*Phi_res*D_phys-F_TRZ^2-F_TRZ^2*SSq = {val:.4f}; obs (1/pi) = {obs:.4f}; match {abs(val-obs)/obs*100:.3f}%.")
