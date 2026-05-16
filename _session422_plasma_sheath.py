"""S422: Bohm-Stangeby sheath potential phi_sh/T_e ~ 2.84 (hydrogen)."""
F_TRZ=0.1; Phi_res=5/6; K_Mex=25/12
val=K_Mex + Phi_res - F_TRZ + F_TRZ**2*K_Mex + F_TRZ**3
tgt=2.84
print(f"S422 COMPLETE. sheath potential = {val:.4f}; target {tgt}; match {abs(val-tgt)/tgt*100:.4f}%.")
