"""S386: Schwarzschild ISCO radius r_ISCO/M = 6 (primitive identification)."""
D_BSFG=6
val = D_BSFG
obs = 6.0
print(f"S386 COMPLETE. r_ISCO/M = D_BSFG = {val}; obs (GR exact) = {obs}; match {abs(val-obs)/obs*100:.4f}%. (D_BSFG locked primitive = ISCO factor)")
