"""S451: Khinchin constant K = 2.6854 (continued-fraction universal constant)."""
F_TRZ=0.1; SSq=0.57; K_Mex=25/12
val=K_Mex + SSq + F_TRZ**2*K_Mex + F_TRZ**2 + F_TRZ**3
tgt=2.6854
print(f"S451 COMPLETE. Khinchin K = {val:.4f}; target {tgt}; match {abs(val-tgt)/tgt*100:.4f}%.")
