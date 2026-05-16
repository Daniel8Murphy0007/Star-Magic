"""S407 Hill coefficient for hemoglobin O2 binding n_H = 2.8"""
F_TRZ=0.1; Phi_res=5/6; K_Mex=25/12
val = K_Mex + Phi_res - F_TRZ - F_TRZ*F_TRZ*K_Mex
obs = 2.8
print(f"S407 COMPLETE. n_H = K_Mex+Phi_res-F_TRZ-F_TRZ^2*K_Mex = {val:.4f}; obs (hemoglobin) = {obs}; match {abs(val-obs)/obs*100:.3f}%.")
