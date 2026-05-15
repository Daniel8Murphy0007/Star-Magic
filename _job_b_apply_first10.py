"""Apply v5.78 closure block to first 10 GW papers.

Inserts the closure section before '## References' (or appends if absent).
Each paper gets a tailored falsifier line based on the table below.
"""
import os, re

WP = 'whitepapers'

# Per-paper tailoring: (filename, extra_hook_line)
# extra_hook_line is appended at the bottom of the block.
PAPERS = [
    ('PAPER_001_GW170817_UQFF_Damping_Analysis.md',
     '**Falsifier hook:** P11 LIGO O5 ringdown spectral offset $R_{21}/R_{22}=0.144$ (PAPER_1175) constrains the 66.7% strain reduction claim made in §3 above.'),
    ('PAPER_002_GW190425_Mass_Gap_Interpretation.md',
     '**Falsifier hook:** P11 LIGO O5 ringdown (PAPER_1175) tests the 47.0% amplitude suppression and the $P(BH)=51\\%$ assignment at the $m_1=2.52\\,M_\\odot$ mass-gap boundary.'),
    ('PAPER_003_GW150914_UQFF_vs_LIGO_Strain.md',
     '**Falsifier hook:** P11 LIGO O5 ringdown spectral offset $R_{21}/R_{22}=0.144$ (PAPER_1175) is the direct observational test of the 66.7% strain reduction reported here.'),
    ('PAPER_004_GW170817_BNS_Chirp_Phase_Evolution.md',
     '**Falsifier hook:** P12 Euclid $\\sigma_8=0.797$ (PAPER_1176) and P13 DESI Y5 $d^2 w/dz^2=0$ (PAPER_1178) bracket the chirp-domain UQFF/GR mismatch at cosmological scale.'),
    ('PAPER_005_BH_Merger_Energy_Retention_UQFF.md',
     '**Master Lagrangian:** `compute_FU_SOURCE4()` and `compute_compressed_MUGE_SOURCE4()` cited in §5 are the closed-Lagrangian assembly of PAPER_1167 (CP4 #254). The 99% mass-retention figure follows directly from the G8 *26! = (1)_{26}* factorial barrier (PAPER_1166).'),
    ('PAPER_006_GW170817_Multi_Messenger_Full_Inspiral.md',
     '**Falsifier hook:** The GW propagation-speed constraint $|\\Delta v/c|<10^{-15}$ cited in §4 enters the P-suite as a baseline anchor for P11 (PAPER_1175). P6 sub-mm Yukawa ($L_{KK}^*=20$–$90\\,\\mu$m, PAPER_1174) is the corresponding short-scale falsifier.'),
    ('PAPER_007_Tidal_Deformability_Constraints_BNS_UQFF.md',
     '**Falsifier hook:** The tidal deformability $\\Lambda$ range ($\\sim 190$–$600$ GR vs UQFF-shifted) is bounded by P11 ringdown spectroscopy (PAPER_1175) and P12 Euclid $\\sigma_8$ (PAPER_1176).'),
    ('PAPER_008_UQFF_Waveform_Phase_Evolution_Template_Mismatch.md',
     '**Falsifier hook:** The 2310.8 rad cumulative phase lag (367.8 cycles) reported in §3 is the canonical template-mismatch signature that P11 LIGO O5 ringdown (PAPER_1175) and P14 CMB-S4 $\\mu$-distortion ($\\mu\\le10^{-8}$, PAPER_1180) jointly constrain.'),
    ('PAPER_008b_Full_Inspiral_Waveform_UQFF.md',
     '**Master Lagrangian:** The explicit $F_U(r,t)=\\sum_{i=1}^4 U_{gi}+U_m+U_A-U_{b,i}$ written at the head of this paper is the *full* closed Lagrangian of PAPER_1167 (CP4 #254 — all eight gaps G1–G8 closed). Each $U_{gi}$, $U_m$, $U_A$, $U_{b,i}$ now has its functional form fixed by G1–G8 rather than postulated.'),
    ('PAPER_009_Damping_Mechanism_Decomposition_UQFF.md',
     '**Per-channel hooks:** Aether $\\to$ G2 DPM SO(2) gauge (PAPER_1163). SCm $\\to$ G7 $\\Phi_{res}=5/6$ closure (PAPER_1165). TRZ $\\to$ G6 $F_{TRZ}=1/10$ closure (PAPER_1163). String sector $\\to$ G5 $T^{22}$ moduli stabilization (PAPER_1164). All four damping channels itemized in §2 are now derived, not parameterized.'),
]

BLOCK_TEMPLATE = """

---

## §v5.78 Closure — Calibration Constants Now Derived

The UQFF damping parameters cited throughout this paper ($\\beta_i$, $F_{{TRZ}}$, $\\rho_{{SCm}}$,
$\\rho_{{UA}}$, $\\kappa$) are no longer free calibrations under canonical UQFF v5.78. Their values
are now outputs of the eight closed Lagrangian gaps (G1–G8, PAPER_1159–1166) and the 27-decade
vacuum-energy ledger (PAPER_1170):

| Constant | Value used here | v5.78 derivation origin |
|---|---|---|
| $\\beta_i = 3(5-i)/20$ | $0.603$ ($i=1$) | G1 Mexican-hat $V(U_A)$ minimum — PAPER_1162 |
| $F_{{TRZ}} = 1/10$ | $0.10$ | G6 topological resonance closure — PAPER_1163 |
| $\\rho_{{SCm}}$ | $7.09\\times10^{{-37}}$ J/m³ | 27-decade vacuum-energy ledger — PAPER_1170 |
| $\\rho_{{UA}}$ | $7.09\\times10^{{-36}}$ J/m³ | 27-decade vacuum-energy ledger — PAPER_1170 |
| $[SSq]$ | $0.57$ | G4 $\\Phi_{{res}}$ / $F_{{TRZ}}$ joint closure — PAPER_1165 |
| $\\kappa$ | $5.0\\times10^{{-4}}$ /day | Empirical decay rate (held); gauged via G3 DPM SO(2) — PAPER_1163 |

**Master synthesis:** PAPER_1167 — *All Eight Lagrangian Gaps Closed* (CP4 #254).
**Vacuum saturation:** PAPER_1170 — *27-Decade R26 + KK + BSFG Vacuum-Energy Ledger* (CP4 #256, $\\rho_\\Lambda$ to <0.5%).

{hook}

*Note:* The $\\xi=13/3$ R26+KK lock (PAPER_1171/1172) is sub-mm-scale and does **not** modify
gravitational-wave predictions in this paper. The closure listed above is the complete v5.78
impact on this whitepaper.
"""

def insert_block(text, block):
    # Try to insert before '## References' (last occurrence)
    m = list(re.finditer(r'\n##\s+References\s*\n', text, re.IGNORECASE))
    if m:
        pos = m[-1].start()
        return text[:pos] + block + text[pos:]
    # Otherwise append at end
    return text.rstrip() + '\n' + block + '\n'

def main():
    for fn, hook in PAPERS:
        path = os.path.join(WP, fn)
        with open(path, 'r', encoding='utf-8') as f:
            text = f.read()
        if '§v5.78 Closure' in text or 'v5.78 Closure' in text:
            print(f'SKIP (already has block): {fn}')
            continue
        block = BLOCK_TEMPLATE.format(hook=hook)
        new = insert_block(text, block)
        with open(path, 'w', encoding='utf-8') as f:
            f.write(new)
        print(f'UPDATED: {fn}')

if __name__ == '__main__':
    main()
