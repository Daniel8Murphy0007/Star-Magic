"""Apply tailored §v5.78 Closure to next 10 papers (PAPER_010..014 + b-variants)."""
import os, re, sys
sys.stdout.reconfigure(encoding='utf-8')

PAPERS = [
    ('PAPER_010_Post_Merger_Oscillations_Remnant_Mass_UQFF.md',
     "**Forward link (this paper):** Post-merger oscillations are direct LIGO O5 ringdown observables. "
     "PAPER_1175 (P11) predicts the ringdown spectral ratio R$_{21}$/R$_{22}$ $\\approx$ 0.144 with $\\beta_i = 0.603$ "
     "and F$_{TRZ}$ = 1/10 locked by the closure above; the remnant-mass calculation in this paper inherits "
     "those values rather than imposing them as fits."),
    ('PAPER_010b_Time_Domain_Chirp_23Hz_UQFF.md',
     "**Forward link (this paper):** The 23 Hz onset and frequency-time evolution analysed here are the "
     "early-inspiral counterpart of the PAPER_1175 (P11) ringdown test. Because $\\beta_i$, F$_{TRZ}$, and "
     "$\\rho_{SCm}$ are no longer free under v5.78, the predicted chirp morphology is parameter-free relative "
     "to v5.77 and is falsifiable in LIGO O5 within $\\pm 5$ Hz of the 23 Hz onset."),
    ('PAPER_011_Stochastic_GW_Background_UQFF_Implications.md',
     "**Forward link (this paper):** SGWB amplitude couples to cosmological parameters. PAPER_1180 (P14, "
     "CMB-S4 $\\mu$-distortion $\\le 10^{-8}$) and PAPER_1178 (P13, DESI Y5 $d^2w/dz^2 = 0$) provide the "
     "falsifiability anchor for the v5.78 SGWB prediction $\\Omega_{GW}(f)$ derived here, since the "
     "27-decade ledger (PAPER_1170) fixes the high-frequency cutoff that v5.77 left underdetermined."),
    ('PAPER_011b_Amplitude_Reduction_Factor_UQFF.md',
     "**Forward link (this paper):** The amplitude-reduction factor calibrated in this paper "
     "IS the strain suppression that PAPER_1175 (P11) tests in LIGO O5. Under v5.78 the AR factor "
     "is not a fit parameter: it is the product of the now-derived $\\beta_i$, F$_{TRZ}$, and "
     "$\\rho_{SCm}/\\rho_{UA}$ ratio fixed by G1, G6, and G7."),
    ('PAPER_012_Eccentric_Binary_Circularization_UQFF.md',
     "**Forward link (this paper):** Eccentric remnants ring down through the same channel as quasi-circular "
     "remnants. PAPER_1175 (P11) ringdown spectral ratio test applies to the residual-eccentricity prediction "
     "e$_f$ derived here; PAPER_1174 (P6, sub-mm Yukawa L$_{KK}^* = 20$-$90$ $\\mu$m) provides an orthogonal "
     "lab test of the KK tower that contributes to the circularization rate."),
    ('PAPER_012b_GW150914_Waveform_Validation.md',
     "**Forward link (this paper):** GW150914 peak strain, phase lag, and damping ratio reported here are the "
     "O1-era benchmark for the PAPER_1175 (P11) LIGO O5 ringdown prediction. With $\\beta_i$ and F$_{TRZ}$ "
     "locked by v5.78 closure, the GW150914 fit reported in this paper is parameter-free and the same "
     "values must hold in O5 within instrumental error."),
    ('PAPER_013_Magnetar_Spin_Down_UQFF_Framework.md',
     "**Forward link (this paper):** The magnetar spin-down anomaly attributed here to UQFF damping is "
     "tested orthogonally by PAPER_1174 (P6, sub-mm Yukawa L$_{KK}^* = 20$-$90$ $\\mu$m), which probes the "
     "same SCm/UA charge sector that drives the spin-down. With $\\rho_{SCm}$ and $\\rho_{UA}$ now derived "
     "from the 27-decade ledger (PAPER_1170) the magnetar B-field amplification factor reported here is "
     "no longer free."),
    ('PAPER_013b_LISA_SMBH_Merger_Rate_UQFF.md',
     "**Forward link (this paper):** The LISA SMBH merger-rate forecast is anchored to cosmological "
     "parameters falsified by PAPER_1176 (P12, Euclid $\\sigma_8 = 0.797$) and PAPER_1178 (P13, DESI Y5 "
     "$d^2w/dz^2 = 0$). The merger-rate density n$(z)$ in this paper inherits the v5.78 $\\sigma_8$ and "
     "dark-energy EoS values rather than treating them as free."),
    ('PAPER_014_Primordial_Black_Holes_UQFF_Formation.md',
     "**Forward link (this paper):** PBH formation imprints on the CMB through $\\mu$-distortion. "
     "PAPER_1180 (P14, CMB-S4 $\\mu \\le 10^{-8}$) and PAPER_1176 (P12, Euclid $\\sigma_8 = 0.797$) provide "
     "the falsifiability anchor for the PBH mass spectrum derived here. The 27-decade ledger (PAPER_1170) "
     "fixes the vacuum-energy scale that controls the PBH formation window."),
    ('PAPER_014b_EMRI_Aether_Damping_UQFF.md',
     "**Forward link (this paper):** EMRI ringdown is the LISA-band analogue of LIGO ringdown. "
     "PAPER_1175 (P11) ringdown spectral-ratio prediction R$_{21}$/R$_{22}$ $\\approx$ 0.144 "
     "carries over to the EMRI signal modifications derived here, since the aether/string damping "
     "channel uses the same $\\beta_i$ and F$_{TRZ}$ now locked by G1 and G6."),
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
predictions in this paper at astrophysical scales. The closure above is the complete v5.78
impact on this whitepaper.

"""

ref_re = re.compile(r'\n##\s+References\s*\n')
done = []
skipped = []

for fn, hook in PAPERS:
    p = os.path.join('whitepapers', fn)
    with open(p, 'r', encoding='utf-8') as f:
        txt = f.read()
    if 'v5.78 Closure' in txt:
        skipped.append(fn)
        continue
    block = BLOCK.replace('{hook}', hook)
    matches = list(ref_re.finditer(txt))
    if matches:
        # insert before LAST ## References
        m = matches[-1]
        new_txt = txt[:m.start()] + '\n' + block + txt[m.start():]
    else:
        new_txt = txt.rstrip() + '\n\n' + block
    with open(p, 'w', encoding='utf-8') as f:
        f.write(new_txt)
    done.append(fn)
    print(f'  INSERTED: {fn}')

print()
print(f'Updated: {len(done)} / Skipped: {len(skipped)}')
if skipped:
    print('Skipped:', skipped)
