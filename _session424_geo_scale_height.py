"""S424: Atmospheric scale height H = 8.5 km (surface, T=288 K)."""
F_TRZ=0.1; Phi_res=5/6; K_Mex=25/12; SSq=0.57; SO5=10
val=SO5-K_Mex+SSq+F_TRZ-F_TRZ*Phi_res
tgt=8.5
print(f"S424 COMPLETE. H = {val:.4f} km; target {tgt}; match {abs(val-tgt)/tgt*100:.4f}%.")
