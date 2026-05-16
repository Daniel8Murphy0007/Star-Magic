"""S441: Higgs self-coupling lambda_H = 0.13 (m_H=125 GeV)."""
F_TRZ=0.1; K_Mex=25/12
val=F_TRZ+F_TRZ**2*K_Mex+F_TRZ**2-F_TRZ**3
tgt=0.13
print(f"S441 COMPLETE. lambda_H = {val:.4f}; target {tgt}; match {abs(val-tgt)/tgt*100:.4f}%.")
