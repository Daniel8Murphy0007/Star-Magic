Dp,FT=4,0.1
v=Dp - FT*Dp + FT
print(f"Ocean avg depth km: {v} vs 3.7 -> {'EXACT' if abs(v-3.7)<1e-9 else 'OFF'}")
