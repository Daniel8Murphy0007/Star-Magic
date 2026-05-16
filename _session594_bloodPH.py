DB,SO,Dp=6,10,4
v=DB + (SO*0.1) + (0.1*Dp)
print(f"Blood pH: {v} vs 7.4 -> {'EXACT' if abs(v-7.4)<1e-9 else 'OFF'}")
