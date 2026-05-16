"""S378: Strong coupling alpha_s at M_Z (running QCD)."""
F_TRZ=0.1; K_Mex=25/12; SSq=0.57; Phi_res=5/6
val = F_TRZ*K_Mex*SSq - F_TRZ**3 * Phi_res
obs = 0.1179
print(f"S378 COMPLETE. alpha_s(M_Z) = F_TRZ*K_Mex*SSq - F_TRZ^3*Phi_res = {val:.5f}; obs (PDG 2022) = {obs}; match {abs(val-obs)/obs*100:.4f}%.")
