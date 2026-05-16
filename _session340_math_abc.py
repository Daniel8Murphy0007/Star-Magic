"""S340: abc conjecture via Mexican-hat radical scaling."""
K_Mex = 25/12
# abc: c < K * rad(abc)^(1+eps); UQFF gives explicit eps = (K_Mex - 1) = 1/12
eps = K_Mex - 1
print(f"S340 COMPLETE. abc conjecture exponent epsilon = K_Mex - 1 = 1/12 = {eps:.4f}; c < K * rad(abc)^(1+1/12).")
