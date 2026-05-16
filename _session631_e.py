from fractions import Fraction as F
F_TRZ=F(1,10);Phi=F(5,6);SSq=F(57,100);K=F(25,12);Dp=4;DB=6;SO=10;N=9;A=60
v=Phi+Phi-F_TRZ**2*DB-F_TRZ**2*SSq
tgt=F(1602,1000)
err=abs(v-tgt)/tgt*100
print(f"Elementary charge e lead: {float(v):.6f} vs {float(tgt):.4f} -> {float(err):.4f}%")
