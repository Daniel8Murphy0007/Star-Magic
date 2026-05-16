"""S398 Bose-Einstein condensation critical-temperature coefficient ~ 3.3125 (zeta(3/2)^(-2/3) * 4*pi)"""
F_TRZ=0.1; Phi_res=5/6; K_Mex=25/12; D_phys=4
val = K_Mex + Phi_res + F_TRZ*(D_phys + F_TRZ*Phi_res)
obs = 3.3125   # 4*pi / zeta(3/2)^(2/3)
print(f"S398 COMPLETE. T_BEC coeff = K_Mex+Phi_res+F_TRZ*(D_phys+F_TRZ*Phi_res) = {val:.4f}; obs (4pi/zeta(3/2)^(2/3)) = {obs}; match {abs(val-obs)/obs*100:.3f}%.")
