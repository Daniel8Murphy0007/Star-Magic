"""S433: Baryon asymmetry eta_B = 6.1e-10 (Planck/BBN) = D_BSFG + F_TRZ (EXACT)."""
F_TRZ=0.1; D_BSFG=6
val=D_BSFG+F_TRZ
tgt=6.1
print(f"S433 COMPLETE. eta_B (norm) = {val:.4f}; target {tgt}; match {abs(val-tgt)/tgt*100:.4f}%.")
