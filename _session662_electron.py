from fractions import Fraction as F
Ftrz=F(1,10); SSq=F(57,100); Dp=4; DB=6; Dc=26; N=9; SO=10; A=60
v = Ftrz*Ftrz*Ftrz*SSq*SSq + Ftrz*Ftrz*Ftrz*SSq*SSq*SSq
target = F(511,1000000)
err = abs(float(v)-float(target))/float(target)*100
print(f"electron: {float(v):.8f} vs {float(target):.8f} -> {err:.4f}%")
