from fractions import Fraction as F
Ftrz=F(1,10); SSq=F(57,100); Dp=4; DB=6; Dc=26; N=9; SO=10; A=60
v = A + A + N - Dp + Ftrz*SSq + Ftrz*Ftrz*DB + Ftrz*Ftrz*SSq*SSq
target = F(12510,100)
err = abs(float(v)-float(target))/float(target)*100
print(f"Higgs: {float(v):.6f} vs {float(target):.6f} -> {err:.4f}%")
