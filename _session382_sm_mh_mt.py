"""S382: Higgs/Top mass ratio m_H/m_t."""
F_TRZ=0.1; K_Mex=25/12; SSq=0.57; beta_i=0.6029
val = beta_i + F_TRZ*K_Mex*SSq
obs = 125.25 / 172.69
print(f"S382 COMPLETE. m_H/m_t = beta_i + F_TRZ*K_Mex*SSq = {val:.4f}; obs (PDG) = {obs:.4f}; match {abs(val-obs)/obs*100:.3f}%.")
