from fractions import Fraction as F
F_TRZ=F(1,10);Phi=F(5,6);SSq=F(57,100);K=F(25,12);Dp=4;DB=6;SO=10;N=9;A=60
v=SO*SSq-F_TRZ**2*Dp+F_TRZ**2
tgt=F(567,100)
err=abs(v-tgt)/tgt*100
print(f"Stefan-Boltzmann lead: {float(v):.6f} vs {float(tgt):.4f} -> {float(err):.4f}%")
