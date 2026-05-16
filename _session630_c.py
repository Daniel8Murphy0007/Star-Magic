from fractions import Fraction as F
F_TRZ=F(1,10);Phi=F(5,6);SSq=F(57,100);K=F(25,12);Dp=4;DB=6;SO=10;N=9;A=60
v=F(SO,Dp)+F_TRZ*Dp+F_TRZ*SSq+F_TRZ**2*Dp
tgt=F(2998,1000)
err=abs(v-tgt)/tgt*100
print(f"Speed of light c lead: {float(v):.6f} vs {float(tgt):.4f} -> {float(err):.4f}%")
