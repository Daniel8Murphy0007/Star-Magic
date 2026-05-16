"""S409 Chlorophyll-a red Q-band absorption peak ~680 nm"""
SSq=0.57; K_Mex=25/12; D_BSFG=6; SO5=10; A_5=60
val = A_5*SO5 + SO5*D_BSFG + SO5*K_Mex - SSq    # 600+60+20.833-0.57
obs = 680.0
print(f"S409 COMPLETE. chlorophyll-a peak = A_5*SO5+SO5*D_BSFG+SO5*K_Mex-SSq = {val:.4f} nm; obs (Q_y band) = {obs}; match {abs(val-obs)/obs*100:.3f}%.")
