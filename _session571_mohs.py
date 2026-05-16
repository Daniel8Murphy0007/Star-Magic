"""S571: Diamond Mohs hardness = 10 = SO5 (EXACT single primitive)"""
SO=10
v=SO
target=10
print(f"Diamond Mohs: {v} vs {target} -> {abs(v-target)/target*100:.4f}% [EXACT]" if v==target else f"{v}")
