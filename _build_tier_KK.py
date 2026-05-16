"""Tier KK builder: generate S684-S693 session scripts from _tier_KK_search.py results."""
from fractions import Fraction

PRIMS = {"F":Fraction(1,10),"Phi":Fraction(5,6),"SSq":Fraction(57,100),"K":Fraction(25,12),"beta":Fraction(6029,10000)}

# Hand-transcribed from _tier_KK_search.py output
ENTRIES = [
    # (ID, name, target, expr_tokens, observed_label)
    (684,"AU_per_1e10_m",          14.96,   ["+K4","-beta","-beta3","-3"],         "1 AU / 10^10 m"),
    (685,"Sun_mass_per_1e29_kg",   19.89,   ["+K4","-2","+3"],                      "M_sun / 10^29 kg"),
    (686,"Earth_orbit_v_per_km_s", 29.78,   ["+K5","-F.K5","-beta","-5"],          "Earth orbital v (km/s)"),
    (687,"Sun_radius_per_1e8_m",   6.96,    ["-F2.beta2","-F2.beta3","+2","+5"],   "R_sun / 10^8 m"),
    (688,"Jupiter_mass_per_1e27_kg",1.898,  ["-F.beta","-F.beta3","-F.beta4","+2"],"M_J / 10^27 kg"),
    (689,"Earth_radius_per_1e6_m", 6.371,   ["+beta2","-2","+3","+5"],             "R_earth / 10^6 m"),
    (690,"Moon_orbital_period_per_day",27.32,["+K5","-F.K5","-3","-5"],            "Moon sidereal period (day)"),
    (691,"sidereal_year_per_100_day",3.65256,["+beta2","+beta3","+F.beta","+3"],   "sidereal year / 100 day"),
    (692,"Mars_orbit_AU",          1.524,   ["-beta2","-beta5","-F.beta2","+2"],   "Mars semi-major axis (AU)"),
    (693,"Mercury_year_per_10_day",8.797,   ["+beta","+beta3","+3","+5"],          "Mercury year / 10 day"),
]

def eval_token(tok):
    sign = Fraction(1,1) if tok[0]=="+" else Fraction(-1,1)
    body = tok[1:]
    # integer literal
    if body.lstrip("-").isdigit(): return sign*Fraction(int(body),1)
    # composite via '.'
    parts = body.split(".")
    prod = Fraction(1,1)
    for p in parts:
        # split primitive name and trailing integer exponent
        i=0
        while i<len(p) and not p[i].isdigit(): i+=1
        nm = p[:i]; ex = int(p[i:]) if i<len(p) else 1
        prod *= PRIMS[nm]**ex
    return sign*prod

TEMPLATE = '''"""S{id} (Tier KK) — {label}.

Locked-primitive closure derived by _tier_KK_search.py.
expr = {expr}
"""
from fractions import Fraction
F=Fraction(1,10); Phi=Fraction(5,6); SSq=Fraction(57,100); K=Fraction(25,12); beta=Fraction(6029,10000)
val = {body}
pred = float(val)
obs  = {obs}
err  = abs(pred-obs)/obs*100
print(f"{name}: {{pred}} vs {{obs}} -> {{err:.4f}}%")
'''

def render_body(tokens):
    parts=[]
    for tok in tokens:
        sign = "+" if tok[0]=="+" else "-"
        body = tok[1:]
        if body.lstrip("-").isdigit():
            parts.append(f"{sign}Fraction({int(body)},1)")
            continue
        segs = body.split(".")
        seg_strs=[]
        for s in segs:
            i=0
            while i<len(s) and not s[i].isdigit(): i+=1
            nm=s[:i]; ex=int(s[i:]) if i<len(s) else 1
            seg_strs.append(f"{nm}**{ex}" if ex!=1 else nm)
        parts.append(f"{sign}({'*'.join(seg_strs)})")
    return "".join(parts)

import pathlib
ROOT = pathlib.Path(__file__).resolve().parent

for sid,name,target,toks,lbl in ENTRIES:
    val = sum((eval_token(t) for t in toks), Fraction(0))
    src = TEMPLATE.format(
        id=sid, label=lbl, expr=" ".join(toks),
        body=render_body(toks),
        obs=target, name=name,
    )
    out = ROOT / f"_session{sid}_KK_{name}.py"
    out.write_text(src, encoding="utf-8")
    pred = float(val); err=abs(pred-target)/target*100
    print(f"S{sid}: {name:32s} {pred:.6g} vs {target}  err={err:.4f}%  -> {out.name}")
