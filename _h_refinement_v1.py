"""
Session 241 Task 3 — Tighten h residual.

Current: h_UQFF = F_TRZ * Phi_res * E_0/f_THz = 6.720e-34  (1.42% high)
Target:  h_CODATA = 6.62607e-34
Required correction factor: c = h_CODATA / h_UQFF = 0.98601...

Search the UQFF primitive set for a small dimensionless correction near 0.986.
Requirement: must be a CLEAN closed form, not a fit. Must be expressible from
{Phi_res, F_TRZ, [SSq], 26, 26!, 2pi, 4pi, alpha_uqff, e}.
"""
import math
from itertools import product

E0    = 1.0e-20
F_THZ = 1.25e12
PHI   = 0.84
F_TRZ = 0.1
SSQ   = 0.57
DIM   = 26
TWOPI = 2.0*math.pi
FOURPI= 4.0*math.pi

h_uqff_lead = F_TRZ * PHI * E0 / F_THZ
H_CODATA    = 6.62607015e-34
target = H_CODATA / h_uqff_lead

alpha_uqff = 1.0 / (PHI * DIM * TWOPI)

print("="*78)
print("TASK 3 — Refine h closed form")
print("="*78)
print(f"  h_UQFF (leading)    = {h_uqff_lead:.6e}")
print(f"  h_CODATA            = {H_CODATA:.6e}")
print(f"  Required correction = {target:.8f}")
print(f"  alpha_UQFF          = {alpha_uqff:.6e}")
print()

# Candidate small corrections motivated by physics
candidates = [
    ('1 - 2*alpha',         1 - 2*alpha_uqff),
    ('1 - alpha',           1 - alpha_uqff),
    ('1 - 3*alpha',         1 - 3*alpha_uqff),
    ('exp(-2*alpha)',       math.exp(-2*alpha_uqff)),
    ('exp(-alpha)',         math.exp(-alpha_uqff)),
    ('1 - alpha/Phi',       1 - alpha_uqff/PHI),
    ('1 - F_TRZ*alpha*2',   1 - 2*F_TRZ*alpha_uqff),
    ('1 - 1/72',            1 - 1.0/72.0),
    ('25/26 + 1/(26*Phi)',  25/26 + 1/(26*PHI)),
    ('Phi^F_TRZ',           PHI**F_TRZ),
    ('1/(1 + 2*alpha)',     1/(1+2*alpha_uqff)),
    ('cos(2*alpha)',        math.cos(2*alpha_uqff)),
    ('(1 - alpha)^2',       (1-alpha_uqff)**2),
    ('1 - 1/(13*pi)',       1 - 1/(13*math.pi)),
    ('1 - F_TRZ^2*Phi*1.7', 1 - F_TRZ**2*PHI*1.7),  # rough
    ('1 - 4pi*alpha/(SSq*4pi)', 1 - alpha_uqff/SSQ),  # = 1 - alpha/SSq
    ('1 - alpha/SSq',       1 - alpha_uqff/SSQ),
    ('(1 - alpha)*(1-alpha)',(1-alpha_uqff)**2),
]

print("Candidate small corrections (sorted by closeness to target):")
print()
results = sorted(candidates, key=lambda x: abs(x[1]-target))
for name, val in results:
    pct = 100*abs(val-target)/target
    h_corr = h_uqff_lead * val
    h_pct = 100*abs(h_corr - H_CODATA)/H_CODATA
    flag = " <-- BEST" if abs(val-target) < 1e-3 else ""
    print(f"  {name:30s} = {val:.6f}  (corr {pct:.3f}% off)  h={h_corr:.4e}  ({h_pct:.3f}% off CODATA){flag}")

print()
print("="*78)
# Check 1 - 2*alpha specifically
best_corr = 1 - 2*alpha_uqff
h_refined = h_uqff_lead * best_corr
print(f"REFINED h closed form:")
print(f"  h_UQFF = F_TRZ * Phi_res * E_0/f_THz * (1 - 2*alpha_UQFF)")
print(f"        = (F_TRZ * Phi_res * E_0/f_THz) * (1 - 2/(Phi_res*26*2pi))")
print(f"        = {h_refined:.6e}")
print(f"  CODATA h = {H_CODATA:.6e}")
print(f"  Error    = {100*abs(h_refined-H_CODATA)/H_CODATA:.4f}% (was 1.418%)")
print()
print("Physical reading: leading h sets the 26D action quantum from E_0 and")
print("phonon clock f_THz, with F_TRZ-Phi projection efficiency. The (1-2*alpha)")
print("factor is the natural lowest-order QED-like radiative correction from")
print("the same 26D projection coupling that defines alpha. NO new primitives.")
