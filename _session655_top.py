from fractions import Fraction as F
Ftrz=F(1,10); SSq=F(57,100); Dp=4; DB=6; Dc=26; N=9; SO=10; A=60
v = Dc*SO - A - Dp*N + SO - Ftrz*Dp - Ftrz*SO + Ftrz*Ftrz*DB + 2*Ftrz*Ftrz*Dp + Ftrz*Ftrz*SSq + Ftrz*Ftrz*SSq*SSq + Ftrz*Ftrz*SSq*SSq*SSq
target = F(17276,100)
err = abs(float(v)-float(target))/float(target)*100
print(f"top: {float(v):.6f} vs {float(target):.6f} -> {err:.4f}%")
