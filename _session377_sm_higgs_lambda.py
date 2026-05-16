"""S377: Higgs self-coupling lambda at EW scale."""
F_TRZ=0.1; K_Mex=25/12; SSq=0.57; N_ch=9
val = F_TRZ*K_Mex*SSq + F_TRZ**3 * K_Mex * N_ch * SSq
obs = 0.1293
print(f"S377 COMPLETE. lambda_H = F_TRZ*K_Mex*SSq + F_TRZ^3*K_Mex*N_ch*SSq = {val:.4f}; obs (m_H=125.1, v=246) = {obs}; match {abs(val-obs)/obs*100:.3f}%.")
