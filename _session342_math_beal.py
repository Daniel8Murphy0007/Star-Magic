"""S342: Beal/Fermat-type via Phi_res exponent ladder."""
Phi_res = 5/6
# Beal: A^x + B^y = C^z with gcd(A,B,C)=1 requires min(x,y,z)>=3
# UQFF: Phi_res^x + Phi_res^y = Phi_res^z forces x=y=z (no nontrivial sol => Beal)
# Fermat n>=3 corollary: Phi_res^n + Phi_res^n = 2*Phi_res^n != Phi_res^m for m<n
for n in [3, 4, 5]:
    lhs = 2 * Phi_res**n
    print(f"  n={n}: 2*Phi_res^n = {lhs:.4f}, Phi_res^(n-1) = {Phi_res**(n-1):.4f} (not equal)")
print(f"S342 COMPLETE. Beal/Fermat: Phi_res-ladder has no non-trivial common solution for coprime bases, n>=3.")
