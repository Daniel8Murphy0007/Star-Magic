Dc,SO,FT=26,10,0.1
v=Dc+SO+FT*SO
print(f"Body temp C: {v} vs 37 -> {'EXACT' if v==37 else 'OFF'}")
