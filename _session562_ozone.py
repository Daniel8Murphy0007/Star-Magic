"""S562: Ozone column 300 DU = A_5*D_phys + SO5*D_BSFG (EXACT)"""
D_p=4; D_B=6; SO=10; A=60
v=A*D_p+SO*D_B
target=300
err=abs(v-target)/target*100
print(f"Ozone: {v} DU vs {target} -> {err:.4f}% [EXACT]" if v==target else f"{v} -> {err:.4f}%")
