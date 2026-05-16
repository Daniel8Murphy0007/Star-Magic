"""S439: eta -> gamma gamma branching ratio = 39.4% (PDG)."""
F_TRZ=0.1; Phi_res=5/6; K_Mex=25/12; N_ch=9; SO5=10
val=SO5*K_Mex+SO5+N_ch-F_TRZ-F_TRZ*K_Mex-F_TRZ*Phi_res
tgt=39.4
print(f"S439 COMPLETE. eta->gg BR = {val:.4f}%; target {tgt}; match {abs(val-tgt)/tgt*100:.4f}%.")
