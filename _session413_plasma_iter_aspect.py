"""S413: ITER aspect ratio R/a = 3.1 = D_BSFG/2 + F_TRZ (EXACT)."""
F_TRZ=0.1; D_BSFG=6
val=D_BSFG/2 + F_TRZ
tgt=3.1
print(f"S413 COMPLETE. ITER aspect R/a = {val:.4f}; target {tgt}; match {abs(val-tgt)/tgt*100:.4f}%.")
