"""S411 Redfield C:N stoichiometric ratio in marine plankton = 106/16 = 6.625"""
F_TRZ=0.1; Phi_res=5/6; D_BSFG=6
val = D_BSFG + Phi_res - F_TRZ - F_TRZ*Phi_res - F_TRZ*F_TRZ - F_TRZ*F_TRZ*F_TRZ
obs = 106.0/16.0
print(f"S411 COMPLETE. Redfield C/N = D_BSFG+Phi_res-F_TRZ-F_TRZ*Phi_res-F_TRZ^2-F_TRZ^3 = {val:.4f}; obs (106/16) = {obs}; match {abs(val-obs)/obs*100:.3f}%.")
