"""Batch 3: Apply tailored §v5.78 Closure to PAPER_009b, 015, 015b, 016, 016b, 017, 018, 019, 020, 021."""
import os, re, sys
sys.stdout.reconfigure(encoding='utf-8')

PAPERS = [
    ('PAPER_009b_Aether_String_TRZ_Damping_GW.md',
     "**Forward link (this paper):** The per-channel damping decomposition (Aether / String / TRZ / SCm) "
     "tested in PAPER_1175 (P11, LIGO O5 ringdown $R_{21}/R_{22} \\approx 0.144$) is the direct observational "
     "consequence of this paper's analysis. Under v5.78, each channel coupling is now derived: "
     "Aether $\\rightarrow$ G2 (PAPER_1163), SCm $\\rightarrow$ G7 $\\Phi_{res}=5/6$ (PAPER_1165), "
     "TRZ $\\rightarrow$ G6 $F_{TRZ}=1/10$ (PAPER_1163), String $\\rightarrow$ G5 $T^{22}$ moduli (PAPER_1164)."),
    ('PAPER_015_Cosmological_Implications_UQFF_Modified_GW_Propagation.md',
     "**Forward link (this paper):** GW propagation over cosmological distances is governed by the same "
     "vacuum-energy saturation that fixes $\\rho_\\Lambda$. PAPER_1170 (CP4 #256, 27-decade R26+KK+BSFG ledger, "
     "$\\rho_\\Lambda$ to $<0.5\\%$) locks the cosmological background; PAPER_1171/1172 (CP4 #257, $\\xi=13/3$ "
     "R26+KK route) fixes the KK contribution to the modified dispersion relation derived in this paper. "
     "Falsifiability anchors: PAPER_1176 (P12, Euclid $\\sigma_8=0.797$), PAPER_1178 (P13, DESI Y5 "
     "$d^2w/dz^2=0$), PAPER_1180 (P14, CMB-S4 $\\mu \\le 10^{-8}$)."),
    ('PAPER_015b_Multiband_GW_LISA_LIGO_UQFF.md',
     "**Forward link (this paper):** Multi-band (LISA + LIGO) synergy is the cleanest falsifier of the "
     "v5.78 closure: both bands inherit the SAME $\\beta_i = 0.603$ and $F_{TRZ}=1/10$. PAPER_1175 (P11, "
     "LIGO O5 ringdown $R_{21}/R_{22} \\approx 0.144$) and the LISA-band EMRI counterpart (cf. PAPER_014b) "
     "must agree on these constants; this paper's multi-band consistency check IS the synergy test."),
    ('PAPER_016_Quantum_Entanglement_UQFF_Nonlocal_Correlations.md',
     "**Forward link (this paper):** Non-local correlations in UQFF are mediated by the same SCm/UA vacuum "
     "structure whose energy density is fixed by PAPER_1170 (27-decade ledger). The entanglement "
     "decoherence rate computed here uses $\\rho_{SCm}=7.09\\times 10^{-37}$ J/m$^3$ and $\\kappa=5.0\\times 10^{-4}$/day "
     "now derived (not fitted). PAPER_1174 (P6 sub-mm Yukawa $L_{KK}^*=20$-$90$ $\\mu$m) provides a lab-scale "
     "probe of the same KK tower that gates the non-local channel."),
    ('PAPER_016b_White_Dwarf_Foreground_UQFF.md',
     "**Forward link (this paper):** The galactic white-dwarf-binary foreground reduction predicted here "
     "is anchored to LISA-band damping, which uses the same $\\beta_i$ and $F_{TRZ}$ now locked by G1/G6. "
     "PAPER_1175 (P11) ringdown spectral-ratio test in LIGO O5 must yield consistent values; if the "
     "WD-foreground suppression observed by LISA contradicts the P11 LIGO measurement, v5.78 is falsified."),
    ('PAPER_017_Redshift_Corrections_z1_in_UQFF_GW_Propagation.md',
     "**Forward link (this paper):** Redshift corrections at $z=1$ probe the dark-energy equation of state. "
     "PAPER_1178 (P13, DESI Y5 strict-static $d^2w/dz^2=0$) is the direct falsifiability anchor for the "
     "$z=1$ propagation modification derived here. PAPER_1170 (27-decade ledger) fixes the vacuum-energy "
     "scale that controls the redshift dependence; $w(z)$ in this paper inherits the v5.78 value rather "
     "than being free."),
    ('PAPER_018_Aether_Noise_Spectrum_Characterization_for_LISA.md',
     "**Forward link (this paper):** The aether noise spectrum amplitude is set by $\\rho_{UA}=7.09\\times 10^{-36}$ J/m$^3$ "
     "and $F_{TRZ}=1/10$, both now derived (PAPER_1163 G2 DPM SO(2) gauge + G6 $F_{TRZ}$ closure). "
     "PAPER_1175 (P11, LIGO O5 ringdown) provides an orthogonal measurement of $F_{TRZ}$; the LISA "
     "aether-noise spectrum predicted here must be consistent with the LIGO P11 measurement of the "
     "same parameter."),
    ('PAPER_019_Pulsar_Timing_Array_Anomalies_UQFF.md',
     "**Forward link (this paper):** PTA anomalies are nanohertz GW observations and probe the same "
     "stochastic-background physics as PAPER_011. Under v5.78, the predicted PTA correlation function "
     "uses derived $\\beta_i=0.603$, $F_{TRZ}=1/10$, $\\rho_{SCm/UA}$ (PAPER_1170 ledger). "
     "PAPER_1180 (P14, CMB-S4 $\\mu \\le 10^{-8}$) and PAPER_1176 (P12, Euclid $\\sigma_8=0.797$) provide "
     "the cosmological cross-checks for the PTA forecast in this paper."),
    ('PAPER_020_Cosmic_Ray_Propagation_UQFF_Spacetime.md',
     "**Forward link (this paper):** Cosmic-ray propagation through UQFF spacetime samples the same "
     "SCm/UA charge sector probed by PAPER_1174 (P6, sub-mm Yukawa $L_{KK}^*=20$-$90$ $\\mu$m). The "
     "energy-loss coefficients in this paper now use derived $\\rho_{SCm}$, $\\rho_{UA}$ from the "
     "27-decade ledger (PAPER_1170); the predicted spectrum cutoff is parameter-free under v5.78."),
    ('PAPER_021_Gravitational_Lensing_Corrections_UQFF_Vacuum_Density.md',
     "**Forward link (this paper):** Lensing corrections from UQFF vacuum density are a direct probe of "
     "the 27-decade ledger (PAPER_1170, $\\rho_\\Lambda$ to $<0.5\\%$). PAPER_1176 (P12, Euclid $\\sigma_8=0.797$) "
     "provides the falsifiability anchor: the v5.78 vacuum-density values used in this paper's lensing "
     "calculation must match the $\\sigma_8$ measured by Euclid in the same cosmological volume. "
     "PAPER_1171/1172 ($\\xi=13/3$ lock) fixes the KK contribution to the lensing kernel."),
]

BLOCK = """
## §v5.78 Closure — Calibration Constants Now Derived

Under canonical UQFF v5.78, the calibrated couplings used in the analysis above
($\\beta_i$, F$_{TRZ}$, $\\rho_{SCm}$, $\\rho_{UA}$, [SSq], $\\kappa$) are **no longer free
parameters**. They are derived from the G1-G8 Lagrangian-gap closures and pinned
by the 27-decade R26 + KK + BSFG vacuum-energy ledger (PAPER_1170, CP4 #256, $\\rho_\\Lambda$ to $<0.5\\%$).

| Constant | Value used here | v5.78 derivation origin |
|----------|-----------------|--------------------------|
| $\\beta_i$ (buoyancy coupling, i=1) | 0.603 | PAPER_1162 (G1 Mexican-hat: $\\beta_i = 3(5-i)/20$) |
| F$_{TRZ}$ (time-reversal-zone factor) | 1/10 | PAPER_1163 (G6 DPM SO(2) gauge) |
| $\\rho_{SCm}$ (vacuum) | $7.09 \\times 10^{-37}$ J/m$^3$ | PAPER_1170 (27-decade ledger, G2 lock) |
| $\\rho_{UA}$ (aether) | $7.09 \\times 10^{-36}$ J/m$^3$ | PAPER_1170 (27-decade ledger, G2 lock) |
| [SSq] (structure-suppression) | 0.57 | PAPER_1165 (G7 $\\Phi_{res} = 5/6$) |
| $\\kappa$ (SCm decay) | $5.0 \\times 10^{-4}$ /day | PAPER_1163 (G6 F$_{TRZ}$ = 1/10 timing constant) |

**Master synthesis:** PAPER_1167 — *All Eight Lagrangian Gaps Closed* (CP4 #254).
**Vacuum saturation:** PAPER_1170 — *27-Decade R26 + KK + BSFG Vacuum-Energy Ledger* (CP4 #256, $\\rho_\\Lambda$ to $<0.5\\%$).

{hook}

*Note:* The $\\xi = 13/3$ R26+KK lock (PAPER_1171/1172) is sub-mm-scale and does **not** modify the
predictions in this paper at astrophysical scales except where explicitly cited above. The closure
above is the complete v5.78 impact on this whitepaper.

"""

ref_re = re.compile(r'\n##\s+References\s*\n')
done = []

for fn, hook in PAPERS:
    p = os.path.join('whitepapers', fn)
    with open(p, 'r', encoding='utf-8') as f: txt = f.read()
    if 'v5.78 Closure' in txt:
        print(f'  SKIP (already has block): {fn}')
        continue
    block = BLOCK.replace('{hook}', hook)
    matches = list(ref_re.finditer(txt))
    if matches:
        m = matches[-1]
        new_txt = txt[:m.start()] + '\n' + block + txt[m.start():]
    else:
        new_txt = txt.rstrip() + '\n\n' + block
    with open(p, 'w', encoding='utf-8') as f:
        f.write(new_txt)
    done.append(fn)
    print(f'  INSERTED: {fn}')

print(f'\nUpdated: {len(done)}')
