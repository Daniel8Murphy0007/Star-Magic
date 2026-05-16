from fractions import Fraction as F
F_TRZ=F(1,10);Phi=F(5,6);SSq=F(57,100);K=F(25,12);Dp=4;DB=6;SO=10;N=9;A=60
v=A*A*Dp*DB + A*SO*N + A*DB*N + SO*N*DB + SO*N*Dp + A*N + Dp + F_TRZ*SO
tgt=96485
err=abs(v-tgt)/tgt*100
print(f"Faraday F: {float(v):.6f} vs {tgt} -> {float(err):.4f}%")
