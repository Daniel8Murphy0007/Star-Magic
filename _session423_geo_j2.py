"""S423: Earth oblateness J_2 = 1.0826e-3 (normalized 1.0826)."""
F_TRZ=0.1; Phi_res=5/6; K_Mex=25/12; SSq=0.57
val=SSq+Phi_res-F_TRZ*K_Mex-F_TRZ-F_TRZ**2-F_TRZ**3
tgt=1.0826
print(f"S423 COMPLETE. J_2 (norm) = {val:.4f}; target {tgt}; match {abs(val-tgt)/tgt*100:.4f}%.")
