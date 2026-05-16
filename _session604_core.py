A,SO,DB,N=60,10,6,9
v=A*SO*DB - SO*SO - DB - N
print(f"Core radius km: {v} vs 3485 -> {'EXACT' if v==3485 else 'OFF'}")
