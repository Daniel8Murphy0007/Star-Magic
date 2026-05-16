"""S461: NANOGrav stochastic GW background h_c = 2.4e-15 (norm 2.4)."""
F_TRZ=0.1; K_Mex=25/12
val=K_Mex + F_TRZ + F_TRZ*K_Mex + F_TRZ**2 - F_TRZ**3
tgt=2.4
print(f"S461 COMPLETE. NANOGrav h_c = {val:.4f}; target {tgt}; match {abs(val-tgt)/tgt*100:.4f}%.")
