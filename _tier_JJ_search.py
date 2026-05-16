"""Brute search for Tier JJ geophysics (Earth + Moon)."""
from fractions import Fraction
F = Fraction(1,10); Phi=Fraction(5,6); SSq=Fraction(57,100); K=Fraction(25,12)
Dp=4; DB=6; Dc=26; N=9; SO=10; A=60; beta=Fraction(6029,10000)

pool=[]; seen=set()
def add(v,l):
    v=Fraction(v)
    if v not in seen: seen.add(v); pool.append((v,l))
for base,bl in [(Dp,'Dp'),(DB,'DB'),(Dc,'Dc'),(N,'N'),(SO,'SO'),(A,'A'),(SSq,'SSq'),(Phi,'Phi'),(K,'K'),(beta,'beta')]:
    for fexp,fl in [(Fraction(1),''),(F,'F'),(F*F,'F2'),(F*F*F,'F3')]:
        for p,pl in [(1,''),(2,'2'),(3,'3'),(4,'4'),(5,'5')]:
            if bl in ('SSq','Phi','beta','K') or p==1:
                v = fexp * (Fraction(base)**p)
                add(v, (fl+'.' if fl else '')+bl+pl)
for k in [1,2,3,4,5]: add(k,str(k))

def search(target, max_terms=6, tol=0.0003):
    target=Fraction(target).limit_denominator(10**8)
    best=None
    n=len(pool)
    def rec(idx,cur,terms,uc):
        nonlocal best
        if uc>0:
            err=abs(float(cur-target)/float(target))
            if best is None or err<best[0]:
                best=(err,float(cur),list(terms))
                if err<tol: return True
        if uc>=max_terms or idx>=n: return False
        if rec(idx+1,cur,terms,uc): return True
        v,l=pool[idx]
        terms.append('+'+l)
        if rec(idx+1,cur+v,terms,uc+1): return True
        terms.pop()
        terms.append('-'+l)
        if rec(idx+1,cur-v,terms,uc+1): return True
        terms.pop()
        return False
    rec(0,Fraction(0),[],0)
    return best

targets=[
 ('Earth_mass_e24kg', 5.972),
 ('Earth_radius_km_e3', 6.371),
 ('Earth_g_ms2', 9.80665),
 ('Earth_escape_kms', 11.186),
 ('Earth_MoI_factor', 0.3307),
 ('Earth_orb_v_kms', 29.78),
 ('Earth_year_days', 365.25),
 ('Moon_dist_e5km', 3.844),
 ('Moon_mass_e22kg', 7.342),
 ('Earth_density_gcm3', 5.514),
]
for name,t in targets:
    b=search(t,max_terms=5,tol=0.0003)
    print(f'{name:22s} target={t}  best={b[1]:.6f}  err={b[0]*100:.4f}%  expr={" ".join(b[2])}')
