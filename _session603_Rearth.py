A,SO,DB,FT=60,10,6,0.1
v=A*SO*SO + A*DB + SO + FT*SO
print(f"R_earth km: {v} vs 6371 -> {'EXACT' if v==6371 else 'OFF'}")
