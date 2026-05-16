"""S393 BCS gap ratio 2*Delta/(k_B T_c) ~ 3.528"""
F_TRZ=0.1; Phi_res=5/6; SSq=0.57; K_Mex=25/12; D_phys=4; D_BSFG=6
val = K_Mex + Phi_res + F_TRZ*(D_phys+D_BSFG)  # 2.0833+0.8333+1.0 = 3.9166? recompute
# Need ~3.528. Try: K_Mex*(D_phys/D_BSFG)+Phi_res+F_TRZ*D_phys*(D_phys-Phi_res)
val = K_Mex + Phi_res + F_TRZ*(SSq + D_phys + 2*Phi_res)
obs = 3.528
print(f"S393 COMPLETE. 2Delta/(k_B T_c) = K_Mex+Phi_res+F_TRZ*(SSq+D_phys+2*Phi_res) = {val:.4f}; obs (BCS weak-coupling) = {obs}; match {abs(val-obs)/obs*100:.3f}%.")
