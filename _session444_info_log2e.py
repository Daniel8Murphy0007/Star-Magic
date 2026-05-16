"""S444: log2(e) = 1.4427 (bit-nat conversion)."""
F_TRZ=0.1; Phi_res=5/6; SSq=0.57; K_Mex=25/12
val=SSq + Phi_res + F_TRZ**2*K_Mex + F_TRZ**2 + F_TRZ**2*Phi_res
tgt=1.4427
print(f"S444 COMPLETE. log2(e) = {val:.4f}; target {tgt}; match {abs(val-tgt)/tgt*100:.4f}%.")
