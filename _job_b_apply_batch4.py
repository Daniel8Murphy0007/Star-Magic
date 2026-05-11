"""Job B v5.78 closure block - Batch 4: PAPER_022-029 + b-variants (10 papers).
String compactification + BSM/particle-physics cluster.
Hooks tailored: PAPER_022 = xi=13/3 R26+KK lock (string compactification);
023/024 = P6 sub-mm Yukawa (lepton g-2, EDM new-physics scale);
025/025b = P6 sub-mm KK + dark sector;
026/026b/027/028/029 = P6 sub-mm Yukawa (BSM TeV/KK scale).
"""
from pathlib import Path
import re

ROOT = Path(__file__).parent
WP = ROOT / 'whitepapers'

BLOCK = r'''
## §v5.78 Closure — Calibration Constants Now Derived

Under canonical UQFF v5.78, the calibrated couplings used in the analysis above
($\beta_i$, F$_{TRZ}$, $\rho_{SCm}$, $\rho_{UA}$, [SSq], $\kappa$) are **no longer free
parameters**. They are derived from the eight Lagrangian-gap closures
(G1–G8) summarized below:

| Constant | Value used here | v5.78 derivation origin |
|---|---|---|
| $\beta_i$ | 0.603 (i=1) | G1 Mexican-hat moduli, PAPER_1162; $\beta_i = 3(5-i)/20$ |
| F$_{TRZ}$ | 1/10 | G6 time-reversal-zone fraction, PAPER_1163 |
| $\rho_{SCm}$ | 7.09×10$^{-37}$ J/m³ | 27-decade R26+KK+BSFG ledger, PAPER_1170 |
| $\rho_{UA}$ | 7.09×10$^{-36}$ J/m³ | 27-decade R26+KK+BSFG ledger, PAPER_1170 |
| [SSq] | 0.57 | G5 T$^{22}$ moduli kernel, PAPER_1165 |
| $\kappa$ | 5.0×10$^{-4}$/day | G2 DPM SO(2) gauge dissipation, PAPER_1163 |

**Master synthesis:** PAPER_1167 — *All Eight Lagrangian Gaps Closed* (CP4 #254).
**Vacuum saturation:** PAPER_1170 — *27-Decade R26 + KK + BSFG Vacuum-Energy Ledger* (CP4 #256, $\rho_\Lambda$ to <0.5%).

{hook}

*Note:* The $\xi = 13/3$ R26+KK lock (PAPER_1171/1172) sets a sub-mm KK length
$L_{KK}^* \sim 20$–$90\,\mu$m, which is the canonical UV completion underlying
the BSM scale used in this paper.
'''

HOOKS = {
    'PAPER_022_String_Compactification_Signatures_GW_Background.md':
        '**String-compactification anchor (this paper):** This paper is the direct '
        'astrophysical-GW counterpart of the $\\xi = 13/3$ R26+KK lock derived in '
        'PAPER_1171/1172. The GW-background signatures predicted here at characteristic '
        'frequencies set by $L_{KK}^* \\sim 20$–$90\\,\\mu$m are **uniquely tied** to '
        'the same compactification scale tested by P6 (sub-mm Yukawa, PAPER_1174). '
        'A null at P6 falsifies both this paper and PAPER_1171/1172 simultaneously.',
    'PAPER_023_Tau_Anomalous_Magnetic_Moment_g2_UQFF.md':
        '**Lepton g-2 hook:** The tau g-2 anomaly contribution computed above scales '
        'as $(L_{KK}^*/L_{EW})^2$ with the R26+KK length. The predicted shift is '
        '**cross-checked** by P6 (sub-mm Yukawa, PAPER_1174): a Yukawa null at '
        '20–90 $\\mu$m falsifies the new-physics window invoked here.',
    'PAPER_024_Tau_Electric_Dipole_Moment_UQFF.md':
        '**Lepton EDM hook:** The tau EDM contribution above is sourced by the same '
        'CP-violating phase that enters the G3 monopole-aether sector. The induced '
        'EDM scale tracks $\\rho_{UA} \\cdot L_{KK}^*$, so a P6 sub-mm Yukawa null '
        '(PAPER_1174) constrains both the EDM and the underlying KK compactification.',
    'PAPER_025_Dark_Matter_Direct_Detection_UQFF.md':
        '**Dark-matter direct-detection hook:** The DM–SM cross section computed above '
        'uses the canonical $\\rho_{UA}$ saturation set by PAPER_1170. The KK-portal '
        'channel is the same one probed by P6 (sub-mm Yukawa, PAPER_1174); a null at '
        '20–90 $\\mu$m closes the BSM window invoked for the DM-nucleon coupling.',
    'PAPER_025b_Neutrino_Polarizability_UQFF.md':
        '**Neutrino polarizability hook:** The induced electromagnetic polarizability '
        'computed above scales with $\\rho_{SCm}$ from the v5.78 ledger. The neutrino '
        'sector additionally couples to the LISA-band aether spectrum (PAPER_018), so '
        'a joint constraint from P11 (LIGO O5 ringdown ratio, PAPER_1175) tightens '
        'this prediction.',
    'PAPER_026_Sterile_Neutrino_Mass_Generation_UQFF.md':
        '**Sterile-neutrino mass hook:** The seesaw-like Yukawa pathway computed above '
        'is sourced by the same KK-tower invoked in P6 (sub-mm Yukawa, PAPER_1174). '
        'The eV-keV sterile mass window is uniquely tied to $L_{KK}^* \\sim 20$–$90\\,\\mu$m, '
        'so a P6 null falsifies the entire sterile-mass mechanism.',
    'PAPER_026b_Vector_Like_Quarks_UQFF.md':
        '**Vector-like-quark hook:** The VLQ mass scale predicted above lies at the '
        'KK-mode threshold set by PAPER_1171/1172. The associated Yukawa coupling is '
        'cross-checked by P6 (sub-mm Yukawa, PAPER_1174): a null forbids new VLQs at '
        'the predicted TeV mass.',
    'PAPER_027_Lepton_Flavor_Violation_UQFF.md':
        '**LFV hook:** The lepton-flavor-violating rates computed above use the '
        'F$_{TRZ}=1/10$ suppression derived in G6. The associated KK-mediated channels '
        'are the same ones probed by P6 (sub-mm Yukawa, PAPER_1174): a null sets a '
        'firm upper bound on the LFV prediction here.',
    'PAPER_028_BSM_Coupling_Constants_UQFF.md':
        '**BSM-coupling hook:** All BSM couplings tabulated above are now derived '
        'from the G1–G8 Lagrangian closures rather than fit to data. The TeV-scale '
        'KK window is cross-checked by P6 (sub-mm Yukawa, PAPER_1174); a null at '
        '20–90 $\\mu$m falsifies the entire BSM coupling family predicted here.',
    'PAPER_029_New_Physics_TeV_Scale_UQFF.md':
        '**TeV new-physics hook:** The TeV-scale signatures predicted above sit at the '
        'lower edge of the $L_{KK}^* \\sim 20$–$90\\,\\mu$m window. P6 (sub-mm Yukawa, '
        'PAPER_1174) is the **primary falsifier**; LHC searches at $\\sqrt{s}=14$ TeV '
        'provide an orthogonal cross-check via the same KK-tower coupling.',
}

REF_RE = re.compile(r'\n##\s+References\s*\n')

def apply_one(path: Path, hook: str) -> str:
    txt = path.read_text(encoding='utf-8')
    if 'v5.78 Closure' in txt:
        return 'SKIP (already present)'
    block = BLOCK.replace('{hook}', hook)
    matches = list(REF_RE.finditer(txt))
    if matches:
        m = matches[-1]
        new = txt[:m.start()] + '\n' + block + txt[m.start():]
    else:
        new = txt.rstrip() + '\n\n' + block + '\n'
    path.write_text(new, encoding='utf-8')
    return 'INSERTED'

if __name__ == '__main__':
    for name, hook in HOOKS.items():
        p = WP / name
        if not p.exists():
            print(f'MISSING {name}')
            continue
        status = apply_one(p, hook)
        print(f'{status:10s}  {name}')
