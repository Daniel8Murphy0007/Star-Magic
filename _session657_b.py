from fractions import Fraction as F
Ftrz=F(1,10); SSq=F(57,100); Dp=4; DB=6; Dc=26; N=9; SO=10; A=60
v = Dp + Ftrz*Dp - Ftrz*SSq - Ftrz*Ftrz*Dc + Ftrz*Ftrz*DB + Ftrz*Ftrz*Dp - Ftrz*Ftrz*SSq*SSq - Ftrz*Ftrz*SSq*SSq*SSq
target = F(418,100)
err = abs(float(v)-float(target))/float(target)*100
print(f"b: {float(v):.6f} vs {float(target):.6f} -> {err:.4f}%")
