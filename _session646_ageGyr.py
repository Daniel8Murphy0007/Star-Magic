from fractions import Fraction as F
Ftrz=F(1,10); Phi=F(5,6); SSq=F(57,100); K=F(25,12)
Dp=4; DB=6; Dc=26; N=9; SO=10; A=60
v = 2*Dp + SO*SSq + Ftrz*SSq + Ftrz*Ftrz*Dp - Ftrz*Ftrz*SSq - Ftrz*Ftrz*SSq*SSq - Ftrz*Ftrz*SSq*SSq*SSq
target = F(1378,100)
err = abs(float(v)-float(target))/float(target)*100
print(f"age_Gyr: {float(v):.6f} vs {float(target):.6f} -> {err:.4f}%")
