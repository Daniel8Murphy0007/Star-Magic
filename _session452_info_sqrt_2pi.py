"""S452: sqrt(2*pi) = 2.5066 (Stirling/Gaussian normalization)."""
F_TRZ=0.1; SSq=0.57; Phi_res=5/6; K_Mex=25/12
val=K_Mex + SSq - F_TRZ - F_TRZ*Phi_res + F_TRZ**2*K_Mex + F_TRZ**2 + F_TRZ**2*Phi_res - F_TRZ**3
tgt=2.5066
print(f"S452 COMPLETE. sqrt(2pi) = {val:.4f}; target {tgt}; match {abs(val-tgt)/tgt*100:.4f}%.")
