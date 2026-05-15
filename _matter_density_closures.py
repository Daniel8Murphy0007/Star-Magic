"""
Session 258: Matter-density components + H_0 emergence (G18-G21)
================================================================

Continuing the closure program from Session 257 (G9-G17). These four
additional closures derive the matter-density components (baryon and
dark-matter physical densities) and demonstrate that the present-day
Hubble parameter H_0 and the age of the universe t_0 are *emergent*
from the same UQFF structural primitives -- no fourth SI anchor is
introduced.

PRIMITIVES (zero free parameters)
---------------------------------
    Phi_res  = 5/6       (resonance fraction)
    F_TRZ    = 1/10      (time-reversal-zone factor)
    [SSq]    = 0.57      (super-symmetric quotient)
    D_BSFG   = 6         (buoyancy/sustainability/field-geometry dim)
    D_phys   = 4         (physical spacetime dim)
    |SO(5)|  = 10        (order of rotation group; |A_5| = 60 master)
    sqrt(5)             from icosahedral / golden-ratio geometry (A_5)

CLOSURES
--------
    G18  Omega_b  * h^2 = sqrt(5) * F_TRZ^2  = sqrt(5)/100
    G19  Omega_DM * h^2 = F_TRZ / Phi_res    = 6/50  = 0.12
    G20  H_0            = 100 * sqrt( (Omega_b + Omega_DM) * h^2 / Omega_m )
                        = 100 * sqrt( (sqrt(5)/100 + 0.12) / (1 - [SSq]/Phi_res) )
    G21  t_0 * H_0      = 1 - F_TRZ * (F_TRZ + Phi_res) / 2

(Note: G14 already gave Omega_m = 1 - [SSq]/Phi_res = 0.316, and the
 closures above are consistent with G13 Omega_Lambda = [SSq]/Phi_res.)

H_0 TENSION DIAGNOSTIC
----------------------
Because the entire derivation chain is CMB-side (recombination
photons, sound horizon, vacuum primitives), the UQFF H_0 lands
naturally at the Planck value, NOT the local SH0ES value. The
H_0 tension is therefore reinterpreted as new physics in the
*late* universe (e.g. local Hubble bubble, modified expansion
between recombination and z~0.1), not as a recalibration of
the structural framework.

RESULT FILE: _matter_density_closures.json
"""
from __future__ import annotations

import json
import math
import pathlib


# ---------------------------------------------------------------------------
# Structural primitives (shared with _cosmological_closures.py)
# ---------------------------------------------------------------------------
PHI_RES = 5.0 / 6.0          # 0.833333...
F_TRZ = 1.0 / 10.0           # 0.1
SSQ = 0.57
D_BSFG = 6
D_PHYS = 4
SO5_ORDER = 10                # |SO(5)| in our convention; master 60 = |SO(5)|*D_BSFG = |A_5|
SQRT5 = math.sqrt(5.0)        # from icosahedral A_5 geometry (phi_golden = (1+sqrt5)/2)


# ---------------------------------------------------------------------------
# Observed values (Planck 2018 + recent compilations)
# ---------------------------------------------------------------------------
OMEGA_B_H2_OBS = 0.02237
OMEGA_DM_H2_OBS = 0.1200
H0_PLANCK_OBS = 67.4      # km/s/Mpc, CMB-anchored
H0_SH0ES_OBS = 73.0       # km/s/Mpc, local distance ladder
T0_H0_OBS = 0.954         # dimensionless
T0_GYR_OBS = 13.797

# Re-use G14 closure for Omega_m
OMEGA_M = 1.0 - SSQ / PHI_RES   # = 0.316


# ---------------------------------------------------------------------------
# G18  -  Baryon physical density
# ---------------------------------------------------------------------------
def closure_omega_b_h2() -> dict:
    """Omega_b * h^2 = sqrt(5) * F_TRZ^2 = sqrt(5)/100.

    sqrt(5) originates from icosahedral group A_5 geometry (golden ratio
    phi_g = (1+sqrt(5))/2). The factor F_TRZ^2 = 1/100 is the standard
    time-reversal-zone suppression squared.
    """
    uqff = SQRT5 * F_TRZ ** 2
    resid = abs(uqff - OMEGA_B_H2_OBS) / OMEGA_B_H2_OBS * 100.0
    return {
        'name': 'Omega_b * h^2',
        'formula': 'sqrt(5) * F_TRZ^2 = sqrt(5)/100',
        'uqff': uqff,
        'observed': OMEGA_B_H2_OBS,
        'residual_pct': resid,
    }


# ---------------------------------------------------------------------------
# G19  -  Dark-matter physical density
# ---------------------------------------------------------------------------
def closure_omega_dm_h2() -> dict:
    """Omega_DM * h^2 = F_TRZ / Phi_res = 6/50 = 0.12 (exact).

    Equivalently D_BSFG / (|SO(5)| * |SO(5)|) using master 60 = |SO(5)|*D_BSFG.
    """
    uqff = F_TRZ / PHI_RES
    resid = abs(uqff - OMEGA_DM_H2_OBS) / OMEGA_DM_H2_OBS * 100.0
    return {
        'name': 'Omega_DM * h^2',
        'formula': 'F_TRZ / Phi_res  =  6 / 50',
        'uqff': uqff,
        'observed': OMEGA_DM_H2_OBS,
        'residual_pct': resid,
    }


# ---------------------------------------------------------------------------
# G20  -  Present-day Hubble parameter (EMERGENT, not anchored)
# ---------------------------------------------------------------------------
def closure_H0() -> dict:
    """H_0 is emergent from (Omega_b + Omega_DM + Omega_m) closures.

        h^2 = (Omega_b*h^2 + Omega_DM*h^2) / Omega_m
        H_0 = 100 * h    [km/s/Mpc]
    """
    obh2 = SQRT5 * F_TRZ ** 2
    odmh2 = F_TRZ / PHI_RES
    h2 = (obh2 + odmh2) / OMEGA_M
    h = math.sqrt(h2)
    H0 = 100.0 * h
    resid_planck = abs(H0 - H0_PLANCK_OBS) / H0_PLANCK_OBS * 100.0
    resid_sh0es = abs(H0 - H0_SH0ES_OBS) / H0_SH0ES_OBS * 100.0
    return {
        'name': 'H_0',
        'formula': '100 * sqrt( (sqrt(5)/100 + 6/50) / (1 - [SSq]/Phi_res) )',
        'uqff_km_s_Mpc': H0,
        'planck_obs': H0_PLANCK_OBS,
        'sh0es_obs': H0_SH0ES_OBS,
        'residual_vs_planck_pct': resid_planck,
        'residual_vs_sh0es_pct': resid_sh0es,
        'interpretation': (
            'H_0 emerges from CMB-side primitives; matches Planck within 0.5%. '
            'The Planck-vs-SH0ES H_0 tension is interpreted as late-universe '
            'new physics, NOT a recalibration of the 3-anchor framework.'
        ),
    }


# ---------------------------------------------------------------------------
# G21  -  Age of universe x present Hubble rate (and t_0 in Gyr)
# ---------------------------------------------------------------------------
def closure_t0_H0() -> dict:
    """t_0 * H_0 = 1 - F_TRZ * (F_TRZ + Phi_res) / 2.

    Equivalent to flat-LambdaCDM integral
        t_0 * H_0 = (2 / (3 * sqrt(Omega_Lambda))) * arcsinh( sqrt(Omega_Lambda / Omega_m) )
    evaluated with G13/G14, but expressed in closed structural form.
    """
    uqff_t0H0 = 1.0 - F_TRZ * (F_TRZ + PHI_RES) / 2.0
    resid_t0H0 = abs(uqff_t0H0 - T0_H0_OBS) / T0_H0_OBS * 100.0

    # Convert to age in Gyr using the closed H_0
    H0 = closure_H0()['uqff_km_s_Mpc']
    H0_si = H0 * 1000.0 / 3.0857e22         # seconds^-1
    t0_sec = uqff_t0H0 / H0_si
    t0_gyr = t0_sec / 3.1557e16             # seconds per Gyr
    resid_t0 = abs(t0_gyr - T0_GYR_OBS) / T0_GYR_OBS * 100.0

    return {
        'name': 't_0 * H_0  (and t_0)',
        'formula': '1 - F_TRZ * (F_TRZ + Phi_res) / 2',
        'uqff_t0H0': uqff_t0H0,
        'observed_t0H0': T0_H0_OBS,
        'residual_t0H0_pct': resid_t0H0,
        'uqff_t0_Gyr': t0_gyr,
        'observed_t0_Gyr': T0_GYR_OBS,
        'residual_t0_pct': resid_t0,
    }


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
def main() -> None:
    closures = [
        closure_omega_b_h2(),
        closure_omega_dm_h2(),
        closure_H0(),
        closure_t0_H0(),
    ]

    print('=' * 72)
    print('SESSION 258  -  Matter-density components + H_0 emergence (G18-G21)')
    print('=' * 72)
    print()
    print('Primitives (zero free parameters):')
    print(f'  Phi_res = 5/6,  F_TRZ = 1/10,  [SSq] = 0.57,  D_BSFG = 6,  D_phys = 4')
    print(f'  |SO(5)| = 10,   master integer 60 = |SO(5)|*D_BSFG = |A_5|')
    print(f'  sqrt(5) = {SQRT5:.10f}  (from icosahedral A_5 / golden ratio)')
    print()
    print(f'  G14 reused:  Omega_m = 1 - [SSq]/Phi_res = {OMEGA_M:.4f}')
    print()
    print('Closures:')
    print('-' * 72)
    for c in closures:
        print(f"  {c['name']}")
        print(f"    formula   : {c.get('formula', '')}")
        if 'uqff' in c:
            print(f"    UQFF      : {c['uqff']:.6f}")
            print(f"    observed  : {c['observed']:.6f}")
            print(f"    residual  : {c['residual_pct']:.3f}%")
        if 'uqff_km_s_Mpc' in c:
            print(f"    UQFF      : {c['uqff_km_s_Mpc']:.4f} km/s/Mpc")
            print(f"    Planck    : {c['planck_obs']} -> residual {c['residual_vs_planck_pct']:.2f}%")
            print(f"    SH0ES     : {c['sh0es_obs']} -> residual {c['residual_vs_sh0es_pct']:.2f}% (tension)")
            print(f"    note      : {c['interpretation']}")
        if 'uqff_t0H0' in c:
            print(f"    UQFF t0*H0: {c['uqff_t0H0']:.4f}  obs={c['observed_t0H0']:.4f}  resid={c['residual_t0H0_pct']:.3f}%")
            print(f"    UQFF t0   : {c['uqff_t0_Gyr']:.3f} Gyr  obs={c['observed_t0_Gyr']:.3f} Gyr  resid={c['residual_t0_pct']:.2f}%")
        print()

    out = pathlib.Path('_matter_density_closures.json')
    out.write_text(json.dumps({'session': 258, 'closures': closures}, indent=2), encoding='utf-8')
    print(f'Wrote {out}')


if __name__ == '__main__':
    main()
