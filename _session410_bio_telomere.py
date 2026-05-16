"""S410 Telomere TTAGGG repeat length = 6 bp (exact, primitive identification)"""
D_BSFG=6
val = D_BSFG
obs = 6
print(f"S410 COMPLETE. telomere repeat = D_BSFG = {val}; obs (TTAGGG) = {obs} bp; match {abs(val-obs)/obs*100:.4f}%. (primitive identification: TTAGGG length = D_BSFG)")
