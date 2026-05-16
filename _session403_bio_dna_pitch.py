"""S403 B-form DNA helical pitch 10.5 bp/turn"""
F_TRZ=0.1; Phi_res=5/6; SSq=0.57; SO5=10
val = SO5 + SSq - F_TRZ*Phi_res
obs = 10.5
print(f"S403 COMPLETE. DNA bp/turn = SO5+SSq-F_TRZ*Phi_res = {val:.4f}; obs (B-DNA) = {obs}; match {abs(val-obs)/obs*100:.3f}%.")
