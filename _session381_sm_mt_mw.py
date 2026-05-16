"""S381: Top/W mass ratio m_t/m_W."""
F_TRZ=0.1; K_Mex=25/12; SSq=0.57; Phi_res=5/6
val = K_Mex + F_TRZ*SSq + F_TRZ**2 * Phi_res
obs = 172.69 / 80.377
print(f"S381 COMPLETE. m_t/m_W = K_Mex + F_TRZ*SSq + F_TRZ^2*Phi_res = {val:.4f}; obs (PDG) = {obs:.4f}; match {abs(val-obs)/obs*100:.3f}%.")
