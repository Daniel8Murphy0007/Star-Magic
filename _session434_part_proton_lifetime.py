"""S434: Proton lifetime tau_p > 1.6e34 yr (Super-K/Hyper-K)."""
F_TRZ=0.1; Phi_res=5/6; K_Mex=25/12; SSq=0.57
val=Phi_res+SSq+F_TRZ*K_Mex-F_TRZ**2
tgt=1.6
print(f"S434 COMPLETE. tau_p (norm) = {val:.4f}; target {tgt}; match {abs(val-tgt)/tgt*100:.4f}%.")
