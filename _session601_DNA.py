SO,Dp,FT=10,4,0.1
v=SO + FT*Dp + FT*FT*SO
print(f"DNA bp/turn: {v} vs 10.5 -> {'EXACT' if abs(v-10.5)<1e-9 else 'OFF'}")
