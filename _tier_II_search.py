"""Brute search for Tier II nuclear binding energies."""
from fractions import Fraction
from itertools import product

F = Fraction(1,10); Phi=Fraction(5,6); SSq=Fraction(57,100); K=Fraction(25,12)
Dp=4; DB=6; Dc=26; N=9; SO=10; A=60; beta=Fraction(6029,10000)

# Build a pool of atomic terms (value, label)
pool = []
def add(v,l): pool.append((Fraction(v),l))
for base,bl in [(Dp,'Dp'),(DB,'DB'),(Dc,'Dc'),(N,'N'),(SO,'SO'),(A,'A'),(SSq,'SSq'),(Phi,'Phi'),(K,'K'),(beta,'beta')]:
    for fexp,fl in [(Fraction(1),''),(F,'F'),(F*F,'F2'),(F*F*F,'F3')]:
        for p,pl in [(1,''),(2,'2'),(3,'3'),(4,'4'),(5,'5')]:
            if bl in ('SSq','Phi','beta','K') or p==1:
                v = fexp * (Fraction(base)**p)
                lbl = (fl+'.' if fl else '') + bl + (pl)
                add(v,lbl)
# also raw small constants
add(1,'1'); add(2,'2'); add(3,'3')

# dedupe
seen=set(); pool2=[]
for v,l in pool:
    if v not in seen:
        seen.add(v); pool2.append((v,l))
pool=pool2

def search(target, max_terms=6, tol=0.001):
    target=Fraction(target).limit_denominator(10**8)
    best=None
    n=len(pool)
    # subset-sum with signs, limit terms
    def rec(idx, cur, terms, used_count):
        nonlocal best
        if used_count>0:
            err=abs(float(cur-target)/float(target))
            if best is None or err<best[0]:
                best=(err, float(cur), list(terms))
                if err<tol: return True
        if used_count>=max_terms or idx>=n: return False
        # skip
        if rec(idx+1,cur,terms,used_count): return True
        v,l=pool[idx]
        terms.append('+'+l)
        if rec(idx+1,cur+v,terms,used_count+1): return True
        terms.pop()
        terms.append('-'+l)
        if rec(idx+1,cur-v,terms,used_count+1): return True
        terms.pop()
        return False
    rec(0,Fraction(0),[],0)
    return best

targets=[
 ('Fe56', 8.7903),
 ('Ni62', 8.7946),
 ('He4', 7.0739),
 ('U235', 7.591),
 ('U238', 7.570),
 ('C12', 7.6802),
 ('O16', 7.9762),
 ('Pb208', 7.8675),
 ('H3', 2.827),
 ('deuteron', 2.2246),
]
for name,t in targets:
    b=search(t,max_terms=5,tol=0.0005)
    print(f'{name:10s} target={t}  best={b[1]:.6f}  err={b[0]*100:.4f}%  expr={" ".join(b[2])}')
