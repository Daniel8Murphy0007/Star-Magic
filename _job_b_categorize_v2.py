#!/usr/bin/env python
"""Job B Phase B1 v2 — content-scan categorizer.

Scans every PAPER_*.md in whitepapers/ for v5.78-relevant content keywords,
and produces a per-paper categorization with explicit reasons.

Buckets (priority C>F>G>D>E>B>A>I>H):
  C  SI constants derivations (alpha, h, c, G, Planck) — filename match
  F  Index / catalog / KB / registry — filename match
  G  Millennium Prize — filename match
  D  KK / extra-dim / sub-mm / Yukawa / compactification — filename or content
  E  Falsifiability / P-suite / DESI / Euclid / CMB-S4 — filename or content
  B  Unified field / Lagrangian / master-equation derivations — filename or content
  A  Cosmology / Lambda / dark energy / vacuum / CMB / BBN — filename or content
  I  LENR / nuclear / Kozima / Holmlid / Pons-Fleischmann — filename, content overrides if ledger/xi present
  H  Specific astrophysical system (catch-all) — content overrides if v5.78 keywords present

For H/I papers a body-content scan promotes them to A/B/D/E if any v5.78
keyword appears in the markdown content (not just filename).
"""

import os, re, csv

WP = 'whitepapers'

# Filename-level rules (high-precision triggers)
FN_C = [r'PAPER_590_', r'PAPER_591_', r'PAPER_592_', r'PAPER_593_', r'PAPER_652_']
FN_F = [r'KnowledgeBase', r'Knowledge_Base', r'_Catalog', r'_Registry', r'_Reference_Table',
        r'_Index', r'Canonical_Body_Reference', r'Species_Index', r'All_8_Lagrangian_Gaps_Closed']
FN_G = [r'Yang_Mills', r'YangMills', r'Navier_Stokes', r'NavierStokes', r'Riemann_Hypothesis',
        r'P_vs_NP', r'PvsNP', r'Hodge_Conjecture', r'Millennium']
FN_D = [r'Compactification', r'Kaluza', r'_KK_', r'Extra_Dimensions', r'Extra_Dim',
        r'Calabi_Yau', r'Type_II[AB]_Superstring', r'String_Theory_.*26D', r'Sub_mm', r'Yukawa']
FN_E = [r'Falsifiab', r'Falsifier', r'_P[0-9]+_', r'DESI', r'Euclid', r'CMB_S4',
        r'LIGO_O5_Ringdown', r'EP1[0-2]_', r'rProcess_Proof']
FN_B = [r'Master_Equation_Derivation', r'Lagrangian_Derivation', r'Lagrangian_Euler_Lagrange',
        r'Wolfram_Field_Unity', r'BSFG_Unification', r'M_Theory_Unification',
        r'UQFF_Phi_Res', r'UQFF_F_TRZ', r'UQFF_26_Factorial', r'UQFF_DPM_SO2',
        r'UQFF_T22_Moduli', r'UQFF_KK_Tower', r'Full_Lagrangian_Unified']
FN_A = [r'Cosmolog', r'Big_Bang', r'BigBang', r'Pre_Big_Bang', r'PreBigBang',
        r'Dark_Energy', r'DarkEnergy', r'Dark_Photon_Manifold',
        r'Vacuum_Energy', r'Vacuum_Density', r'VacuumEnergy', r'VacuumDensity',
        r'Lambda_CDM', r'LambdaCDM', r'LCDM', r'_Lambda_', r'Cosmological_Constant',
        r'Reionization', r'CMB_', r'BBN', r'Inflation', r'Inflationary',
        r'Universe_Diameter', r'UniverseDiameter', r'Observable_Universe',
        r'Cosmogenesis', r'Quintessence',
        r'PAPER_1156_', r'PAPER_1170_', r'PAPER_1119_']
FN_I = [r'_LENR', r'LENR_', r'Kozima', r'Holmlid', r'Pons[_-]?Fleischmann',
        r'Widom_Larsen', r'_WL_', r'Cold_Fusion', r'Muon_Catalyzed']

# Content-scan keywords — HIGH-SPECIFICITY v5.78-only tokens.
# (R26, KK tower, F_TRZ etc. are core UQFF vocab; they appear in hundreds of papers
#  and are NOT v5.78 markers. Only the explicit closure/ledger/xi=13/3 tokens are.)
KW_LEDGER   = [r'27-decade', r'27 decade', r'vacuum[- ]energy ledger',
               r'PAPER_1170', r'27\s*orders of magnitude.{0,40}vacuum']
KW_XI       = [r'\bxi\s*=\s*13/3\b', r'(?<!\d)13/3(?!\d)\s*(lock|locked|saturation|saturates|closure)',
               r'PAPER_1171', r'PAPER_1172',
               r'Gauss-?Bonnet.{0,30}R26', r'R26.{0,30}Gauss-?Bonnet']
KW_LAG_GAP  = [r'G1[-\s]?G8', r'\bG[1-8]\s+gap\b', r'PAPER_116[0-7]',
               r'eight Lagrangian gaps', r'closed Lagrangian', r'Lagrangian gap closure',
               r'Mexican-?hat.{0,30}V\(UA\)', r'beta_?i\s*=\s*3\(5']
KW_SI       = [r'three-?anchor', r'PAPER_590\b', r'PAPER_591\b', r'PAPER_592\b',
               r'PAPER_593\b', r'three constants closure',
               r'alpha\s*=\s*1\s*/\s*\(?\s*26\s*[*xX·]\s*2\s*pi']
KW_FALS     = [r'\bP[6-9]\b\s*(predict|falsif|sub-?mm)', r'\bP1[0-4]\b\s*(predict|falsif)',
               r'PAPER_117[4-9]', r'PAPER_1180',
               r'P1-?P14', r'P6-?P14', r'P6[\-\u2013]P10',
               r'DESI\s*Y5', r'CMB-?S4\s*mu', r'LIGO\s*O5\s*ringdown',
               r'Euclid\s*sigma_?8\s*=']
KW_KK       = [r'KK regulator', r'sub-?mm Yukawa\b',
               r'L_?KK\s*\*', r'PAPER_1173\b', r'L_?KK\s*=\s*20[-\u2013]90']

CONTENT_BUCKETS = [
    ('A', KW_LEDGER, 'content: 27-decade vacuum ledger / Lambda closure'),
    ('B', KW_LAG_GAP, 'content: G1-G8 closed Lagrangian'),
    ('D', KW_KK, 'content: KK tower / sub-mm Yukawa'),
    ('E', KW_FALS, 'content: P-suite falsifiability'),
    ('C', KW_SI, 'content: three-anchor SI closure'),
    # xi=13/3 alone leans D (R26+KK lock paper), or A if also ledger
    ('D', KW_XI, 'content: xi=13/3 R26+KK lock'),
]

ACTIONS = {
    'C': 'Three-anchor SI closure banner (Sess 237-241); STRUCTURAL tier',
    'A': 'Add closing section: v5.78 27-decade vacuum ledger (PAPER_1170) + xi=13/3',
    'B': 'Add closed-Lagrangian cross-ref (PAPER_1159-1167 G1-G8) + CP4 #254',
    'D': 'Forward-pointer to PAPER_1171/1172 (xi=13/3 R26+KK lock) + PAPER_1173',
    'E': 'Forward-pointer to PAPER_1174 P1-P14 suite + PAPER_1177-1180',
    'F': 'Update tables; add P6-P14, CP4 #254-#264, ref PAPER_1167/1170/1174',
    'G': 'Verify v5.78 closed-Lagrangian/ledger cross-ref is present (Sess 225 baseline)',
    'I': 'NO UPDATE unless ledger/xi/KK/P-suite enters body (per Job B doc)',
    'H': 'NO UPDATE (specific system / framework application)',
}

def match_any(patterns, text):
    for p in patterns:
        if re.search(p, text, re.IGNORECASE):
            return p
    return None

def classify_filename(fn):
    for label, pats in [('C',FN_C),('F',FN_F),('G',FN_G),('D',FN_D),
                        ('E',FN_E),('B',FN_B),('A',FN_A),('I',FN_I)]:
        m = match_any(pats, fn)
        if m:
            return label, f'filename: {m}'
    return 'H', 'filename: no v5.78 trigger -> default H'

def classify_content(text):
    for label, kws, reason in CONTENT_BUCKETS:
        m = match_any(kws, text)
        if m:
            return label, f'{reason} ({m})'
    return None, None

def main():
    rows = []
    for fn in sorted(os.listdir(WP)):
        if not fn.startswith('PAPER_') or not fn.endswith('.md'):
            continue
        path = os.path.join(WP, fn)
        try:
            with open(path, 'r', encoding='utf-8', errors='replace') as f:
                text = f.read()
        except Exception as e:
            text = ''
        fn_bucket, fn_reason = classify_filename(fn)
        final_bucket = fn_bucket
        final_reason = fn_reason
        # Content override only for H or I (default catch-alls)
        if fn_bucket in ('H', 'I'):
            c_bucket, c_reason = classify_content(text)
            if c_bucket:
                final_bucket = c_bucket
                final_reason = f'PROMOTED from {fn_bucket} -> {c_bucket} ({c_reason})'
        m = re.match(r'PAPER_(\d+)([a-z]?)', fn)
        pid = m.group(1) if m else '0000'
        variant = m.group(2) if m else ''
        rows.append({
            'paper_id': pid,
            'variant': variant,
            'filename': fn,
            'bucket': final_bucket,
            'reason': final_reason,
            'action': ACTIONS[final_bucket],
        })
    # Sort numerically
    rows.sort(key=lambda r: (int(r['paper_id']), r['variant']))
    out = '_job_b_categorization_v2.csv'
    with open(out, 'w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['paper_id','variant','filename','bucket','reason','action'])
        w.writeheader()
        w.writerows(rows)
    # Counts
    from collections import Counter
    cnt = Counter(r['bucket'] for r in rows)
    print(f'Scanned: {len(rows)} papers')
    print('Bucket counts:')
    for b in ['C','A','B','D','E','F','G','I','H']:
        print(f'  {b}: {cnt.get(b,0)}')
    print(f'NEEDS UPDATE (C+A+B+D+E+F+G): {sum(cnt.get(b,0) for b in "CABDEFG")}')
    print(f'PROMOTED via content scan: {sum(1 for r in rows if r["reason"].startswith("PROMOTED"))}')
    print(f'Wrote {out}')

if __name__ == '__main__':
    main()
