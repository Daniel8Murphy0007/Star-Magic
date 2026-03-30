"""Add Session 153 Alders/Olbers classes (PAPER_564-566) and
   Session 153b gap-fill CP4 classes (PAPER_567-572) to CondensedPhysics4.py.
   Run from repo root. Uses read/write in utf-8 (no BOM).
"""
import pathlib, sys

path = pathlib.Path('CondensedPhysics4.py')
content = path.read_text(encoding='utf-8')

# ------------------------------------------------------------------ #
# Insert block 1: constants + 3 Olbers classes, before __all__
# ------------------------------------------------------------------ #
ANCHOR_BEFORE_ALL = "            'session': self.SESSION, 'papers': [self.PAPER],\n        }\n\n\n__all__ = ["

NEW_CLASSES = r"""            'session': self.SESSION, 'papers': [self.PAPER],
        }


# ─────────────────────────────────────────────────────────────────────────────
# Session 153 constants — Alders/Olbers Paradox Resolution
# ─────────────────────────────────────────────────────────────────────────────
_S153_SSq       = 0.507
_S153_H0        = 2.268e-18          # 70 km/s/Mpc [s^-1]
_S153_c         = 2.998e8            # speed of light [m/s]
_S153_rH        = 4.4e26             # c/H0 [m]
_S153_nstar_m3  = 1.0e9 / (3.086e22 ** 3)   # 1e9 Mpc^-3 -> m^-3
_S153_Lstar     = 3.828e26           # 1 L_sun [W]
_S153_rho_UA    = 7.09e-36           # UA' vacuum [J/m^3]
_S153_rho_SCm   = 7.09e-37           # SCm vacuum [J/m^3]
_S153_kappa_DPM = 5.0e-4
_S153_N         = 26
_S153_phi       = (1.0 + 5.0 ** 0.5) / 2.0
_S153_omega_CW  = 2.0 * 3.141592653589793 * 1.2e10   # PAPER_516
_S153_ETA       = 1.0e-22            # BSFG aether eta
_S153_MS        = 1.989e30           # solar mass [kg]


class AldersOlbersParadoxDPMShellFluxCalculator(_CP4Calculator):
    """CP4 #158 -- PAPER_564: Alders/Olbers Paradox via DPM 26-Shell Radiance Cascade.

    Classical problem: shell-by-shell sky brightness diverges (B_sky -> inf).
    UQFF Resolution: [SSq]-damped DPM radiance per shell.

    Shell n flux: B_n = n_star*L_star*dr_n / (4*pi*c*(1+z_n)^4) * R_{Ug1,n}
    R_{Ug1,n} = F*(1+M_sf)*exp(-[SSq]*n/26)          [PAPER_427]
    z_n = H0*r_n/c,  dr_n = r_H/26
    Total: B_sky = sum_{n=1}^{26} B_n  CONVERGES via [SSq]-geometric damping.
    ShellEnergy^(n) = kappa*(DPM_n-DPM_s)/r_n^26 * omega_CW  [PAPER_516]
    ProtoH: B_total = B_sky + DPM_react*Prob_order*|t_neg|   [PAPER_519]

    Builds on: DPMLayeredShellEnergyRadianceCalculator (#111, PAPER_516)
               TwentySixDResonanceLayerAmplitudeFrequencyCalculator (#81, PAPER_427)
               ShellRadiancePrototypeEquationCalculator (#114, PAPER_519)
    Source: Session 153 -- Alders/Olbers codebase resolution  PAPER_564
    """
    SESSION = 153
    PAPER   = 'PAPER_564'

    def compute(self, dataset: dict = None) -> dict:
        import math
        if dataset is None:
            dataset = {}
        n_star   = dataset.get('n_star',   _S153_nstar_m3)
        L_star   = dataset.get('L_star',   _S153_Lstar)
        H0       = dataset.get('H0',       _S153_H0)
        c        = dataset.get('c',        _S153_c)
        r_H      = dataset.get('r_H',      _S153_rH)
        SSq      = dataset.get('SSq',      _S153_SSq)
        kappa    = dataset.get('kappa',    _S153_kappa_DPM)
        DPM_n    = dataset.get('DPM_n',    1.0)
        DPM_s    = dataset.get('DPM_s',    0.85)
        omega_CW = dataset.get('omega_CW', _S153_omega_CW)
        F_amp    = dataset.get('F_amp',    1.0)
        M_sf     = dataset.get('M_sf',     0.1)
        t_neg    = dataset.get('t_neg',   -1.0)
        N  = _S153_N
        dr = r_H / N
        shells, B_sky = [], 0.0
        for i in range(1, N + 1):
            r_n       = i * dr
            z_n       = H0 * r_n / c
            dimming   = (1.0 + z_n) ** 4
            R_n       = F_amp * (1.0 + M_sf) * math.exp(-SSq * i / N)
            B_n       = (n_star * L_star * dr / (4.0 * math.pi * c)) / dimming * R_n
            shell_e_n = kappa * (DPM_n - DPM_s) / (r_n ** 26 + 1e-300) * omega_CW
            shells.append({'n': i, 'r_n_m': r_n, 'z_n': z_n,
                           'R_Ug1_n': R_n, 'B_n_W_m2_sr': B_n,
                           'shell_energy': shell_e_n})
            B_sky += B_n
        DPM_react  = kappa * (DPM_n - DPM_s) / r_H
        Prob_order = math.exp(-1.0 / 9.0) / (1.0 + abs(t_neg))
        B_DPM_corr = DPM_react * Prob_order * abs(t_neg)
        B_classical = n_star * L_star * r_H / (4.0 * math.pi * c)
        conv_ratio  = B_sky / B_classical if B_classical > 0 else float('nan')
        return {
            'paper': self.PAPER,
            'B_sky_total_W_m2_sr': B_sky,
            'B_classical_divergent': B_classical,
            'convergence_ratio': conv_ratio,
            'B_DPM_correction': B_DPM_corr,
            'N_shells': N,
            'SSq_damping_factor': SSq,
            'r_H_m': r_H, 'dr_m': dr,
            'shells_preview': shells[:4],
            'paradox_resolved': B_sky < B_classical,
            'primary_equations': [
                'B_n = n_star*L_star*dr / (4*pi*c*(1+z_n)^4) * R_{Ug1,n}',
                'R_{Ug1,n} = F*(1+M_sf)*exp(-[SSq]*n/26)  [PAPER_427]',
                'B_sky = sum_{n=1}^{26} B_n  CONVERGES',
                'z_n = H0*r_n/c  (Hubble redshift per shell)',
                'ShellEnergy^(n) = kappa*(DPM_n-DPM_s)/r_n^26 * omega_CW  [PAPER_516]',
            ],
            'available_equations': [
                'B_classical = n*L*r_H/(4*pi*c)  diverges',
                'Convergence ratio = B_sky_UQFF / B_classical',
                'DPM_react = kappa*(DPM_n-DPM_s)/r_H',
                'Prob_order = exp(-S_{26D}/v_init)/Partition_{9D}  [PAPER_519]',
                'B_total = B_sky + DPM_react*Prob_order*|t_neg|  [ProtoH]',
            ],
            'simulation_set': [
                'Sweep [SSq] in [0.3,0.9] -> B_sky convergence envelope',
                'Vary N 1->26 -> monotone convergence of B_sky',
                'Sweep H0 in [67.4,73.0] km/s/Mpc (Hubble tension)',
                'Plot B_n vs n: exponential [SSq]-decay curve',
            ],
            'session': self.SESSION, 'papers': [self.PAPER],
        }


class AldersOlbersVDSNumberSystemResolutionCalculator(_CP4Calculator):
    """CP4 #159 -- PAPER_565: Alders/Olbers Paradox via VDS/DVP/BH Number Systems.

    VDS (Vacuum Density Series) photon flux bound [PAPER_429]:
      B_sky^VDS = sum_{k=1}^{26} (1/k^26)*[SSq]^k * (n*L*dr_k/(4*pi*c))
      Formal bound: B_sky <= n*L*r_H/(4*pi*c) * Li_{26}([SSq])
      Li_{26}([SSq]=0.507) ~ 0.507  (49.3% suppression)

    DVP (Dipole Vortex Primes) vortex scattering [PAPER_429]:
      A(p) proportional to [SSq]^{pi(p)}/p^26  for primes p>26
      p_special=113 (H proto-shell anchor)

    BH (Buoyancy Harmonics) absorption [PAPER_429]:
      H_m = sum_{k=1}^m (1/k)*f_Ub
      U_g2 = sum_m H_m*(1-exp(-[SSq]*m))*cos(omega_Ug2*t_n)

    Unification [PAPER_535]: Z = Li_{26}([SSq]) -- VDS+DVP+BH unified

    Source: Session 153 -- Alders/Olbers codebase resolution  PAPER_565
    """
    SESSION = 153
    PAPER   = 'PAPER_565'

    def compute(self, dataset: dict = None) -> dict:
        import math
        if dataset is None:
            dataset = {}
        SSq       = dataset.get('SSq',       _S153_SSq)
        n_star    = dataset.get('n_star',    _S153_nstar_m3)
        L_star    = dataset.get('L_star',    _S153_Lstar)
        c         = dataset.get('c',         _S153_c)
        r_H       = dataset.get('r_H',       _S153_rH)
        f_Ub      = dataset.get('f_Ub',      1.0e-3)
        omega_Ug2 = dataset.get('omega_Ug2', 1.0e12)
        t_n       = dataset.get('t_n',       0.0)
        N  = _S153_N
        dr = r_H / N
        # Li_26([SSq])
        Li26 = 0.0
        for k in range(1, 201):
            term = (SSq ** k) / (k ** 26)
            Li26 += term
            if term < 1.0e-30:
                break
        B_VDS = sum(
            (SSq ** k) / (k ** 26) * n_star * L_star * dr / (4.0 * math.pi * c)
            for k in range(1, N + 1)
        )
        B_classical  = n_star * L_star * r_H / (4.0 * math.pi * c)
        B_bound_Li26 = B_classical * Li26

        def _is_prime(n):
            if n < 2: return False
            if n == 2: return True
            if n % 2 == 0: return False
            return all(n % i for i in range(3, int(n**0.5)+1, 2))

        primes_gt26 = [p for p in range(27, 200) if _is_prime(p)]
        pi_count    = len(primes_gt26)
        pi_113      = sum(1 for p in range(2, 114) if _is_prime(p))
        A_DVP_total = sum(
            (SSq ** (i + 1)) / (p ** 26)
            for i, p in enumerate(primes_gt26)
        )
        ell_DVP = r_H / max(pi_count, 1)

        two_pi = 2.0 * math.pi
        U_g2_total = 0.0
        for m in range(1, N + 1):
            H_m       = sum(f_Ub / k for k in range(1, m + 1))
            U_g2_total += H_m * (1.0 - math.exp(-SSq * m)) * math.cos(omega_Ug2 * t_n * two_pi)

        SSq_dyn = math.log(_S153_rho_SCm / _S153_rho_UA) * 13.0 * math.exp(-(math.pi - t_n))

        return {
            'paper': self.PAPER,
            'Li26_SSq': Li26,
            'B_classical_W_m2_sr': B_classical,
            'B_bound_Li26_W_m2_sr': B_bound_Li26,
            'B_VDS_W_m2_sr': B_VDS,
            'VDS_suppression_fraction': 1.0 - (B_VDS / B_classical) if B_classical > 0 else float('nan'),
            'Z_unification': Li26,
            'DVP_primes_gt26_first10': primes_gt26[:10],
            'DVP_pi_count': pi_count,
            'DVP_A_total': A_DVP_total,
            'DVP_pi_113': pi_113,
            'DVP_ell_mfp_m': ell_DVP,
            'BH_U_g2_total': U_g2_total,
            'SSq_dynamic_n13': SSq_dyn,
            'paradox_resolved': B_VDS < B_classical,
            'primary_equations': [
                'Z = Li_{26}([SSq]) = sum_{k=1}^inf [SSq]^k/k^26  [PAPER_535]',
                'B_sky <= n*L*r_H/(4*pi*c) * Li_{26}([SSq])       [VDS bound]',
                'B^VDS = sum_{k=1}^{26} [SSq]^k/k^26 * n*L*dr_k/(4*pi*c)',
                'A(p)_DVP proportional to [SSq]^{pi(p)}/p^26  p>26 p_special=113',
                'U_g2 = sum_m H_m*(1-exp(-[SSq]*m))*cos(omega*t_n)',
                '[SSq]_dyn(n,t)=log(rho_SCm/rho_UA)*n*exp(-(pi-t_n))  [PAPER_429]',
            ],
            'available_equations': [
                'Mean free path: ell_DVP = r_H / pi_count(DVP primes 27..200)',
                'H_m = sum_{k=1}^m (1/k)*f_Ub  (BH harmonic)',
                'VDS convergence condition: |[SSq]|<1  satisfied: [SSq]=0.507',
                'B_sky/B_classical = Li_{26}([SSq]) (fractional sky brightness)',
            ],
            'simulation_set': [
                'Vary [SSq] in [0.3,1.0] -> Li_{26} diverges at |[SSq]|->1',
                'DVP primes p in (26,500) -> scattering amplitude spectrum A(p)',
                'BH harmonic absorption U_g2 vs m at THz boundary',
                'Dynamic [SSq](n,t): sweep t_n in [-pi, pi] at n=13',
            ],
            'session': self.SESSION, 'papers': [self.PAPER],
        }


class AldersOlbersBSFGMetricGapAnalysisCalculator(_CP4Calculator):
    """CP4 #160 -- PAPER_566: Alders/Olbers Paradox via BSFG Aether Metric + Gap Analysis.

    BSFG Resolution [PAPER_554]:
      A_uv = g_uv + eta*T_s00*cos(pi*t_n)*delta_uv
      R^r_0r0 = 6*eta*C_num/r^5
      Gamma_BSFG = eta*|R_scalar_avg|/c^4
      E(r) = E_0*exp(-Gamma_BSFG*r)  [photon energy extinction]
      B_n^BSFG = n*L*dr/(4*pi*c) * exp(-Gamma*r_n) * [SSq]^n
      B_sky^BSFG = sum_{n=1}^{26} B_n^BSFG  [VDS x BSFG double suppression]

    Gap Analysis (6 present / 6 missing):
      PRESENT: Hubble finite-age, (1+z)^4, DPM 26-shell, [SSq] damping,
               VDS+Li_{26}, BSFG aether geodesic
      MISSING: n_star(z) SFR evolution, kappa_lambda opacity,
               B_sky_obs benchmark, DVP photon-photon scatter,
               t_neg timing, W/m^2/sr calibration

    Builds on: BSFGRiemannCurvatureAetherMetricCalculator (#149, PAPER_554)
               AldersOlbersParadoxDPMShellFluxCalculator (#158, PAPER_564)
               AldersOlbersVDSNumberSystemResolutionCalculator (#159, PAPER_565)
    Source: Session 153 -- Alders/Olbers codebase resolution  PAPER_566
    """
    SESSION = 153
    PAPER   = 'PAPER_566'

    def compute(self, dataset: dict = None) -> dict:
        import math
        if dataset is None:
            dataset = {}
        SSq    = dataset.get('SSq',    _S153_SSq)
        n_star = dataset.get('n_star', _S153_nstar_m3)
        L_star = dataset.get('L_star', _S153_Lstar)
        c      = dataset.get('c',      _S153_c)
        r_H    = dataset.get('r_H',    _S153_rH)
        eta    = dataset.get('eta',    _S153_ETA)
        Ms     = dataset.get('Ms',     _S153_MS)
        Ls     = dataset.get('Ls',     3.828e26)
        t_n    = dataset.get('t_n',    0.0)

        C_num        = (Ms * c**2 + Ls / c**2) / ((4.0/3.0) * math.pi)
        R_scalar_avg = 6.0 * eta * C_num / (r_H ** 5)
        Gamma_BSFG   = eta * abs(R_scalar_avg) / (c ** 4)

        N  = _S153_N
        dr = r_H / N
        B_sky_BSFG = 0.0
        shells_bsfg = []
        for n in range(1, N + 1):
            r_n     = n * dr
            E_ratio = math.exp(-Gamma_BSFG * r_n)
            vds_fac = SSq ** n
            B_n     = n_star * L_star * dr / (4.0 * math.pi * c) * E_ratio * vds_fac
            shells_bsfg.append({'n': n, 'r_n_m': r_n,
                                 'BSFG_E_ratio': E_ratio,
                                 'VDS_factor': vds_fac, 'B_n_BSFG': B_n})
            B_sky_BSFG += B_n
        B_classical = n_star * L_star * r_H / (4.0 * math.pi * c)

        gap_analysis = {
            'PRESENT': {
                'Hubble_horizon':    'RedshiftDependentHubbleCalculator CP2:15330',
                '1pz4_dimming':      'H(t,z)=H0*sqrt(0.3*(1+z)^3+0.7) CP2',
                'DPM_26shell':       'DPMLayeredShellEnergyRadianceCalculator CP4 #111 PAPER_516',
                'SSq_damping':       'TwentySixDResonanceLayerAmplitudeFrequencyCalculator CP4 #81 PAPER_427',
                'VDS_Li26':          'ThreeNewNumberSystemsVacuumDipoleBuoyancyCalculator+VDSDVPBHNumberSystemsCatalogue PAPER_429+535',
                'BSFG_geodesic':     'BSFGRiemannCurvatureAetherMetricCalculator CP4 #149 PAPER_554',
                'Olbers_DPM':        'AldersOlbersParadoxDPMShellFluxCalculator CP4 #158 PAPER_564',
                'Olbers_VDS':        'AldersOlbersVDSNumberSystemResolutionCalculator CP4 #159 PAPER_565',
            },
            'MISSING': {
                'n_star_z':       'n_star(z) SFR Madau-Dickinson stellar density evolution -- PAPER_567',
                'kappa_lambda':   'kappa_lambda(lambda) wavelength-dependent opacity -- PAPER_568',
                'B_obs_bench':    'B_sky_obs=3.1e-6 W/m^2/sr CMB/EBL benchmark validation -- PAPER_569',
                'DVP_scatter':    'Photon-photon DVP prime vortex encompassment scattering -- PAPER_570',
                't_neg_timing':   't_neg photon arrival timing DPM delay formalisation -- PAPER_571',
                'unit_calib':     'Shell radiance calibrated to observable W/m^2/sr -- PAPER_572',
            },
            'completeness': '6/6 primary resolutions present; 6 extensions missing (PAPER_567-572)',
        }

        return {
            'paper': self.PAPER,
            'Gamma_BSFG_per_m': Gamma_BSFG,
            'R_scalar_avg': R_scalar_avg,
            'C_num': C_num,
            'B_sky_BSFG_W_m2_sr': B_sky_BSFG,
            'B_classical': B_classical,
            'BSFG_suppression': B_sky_BSFG / B_classical if B_classical > 0 else float('nan'),
            'shells_preview': shells_bsfg[:4],
            'gap_analysis': gap_analysis,
            'paradox_resolved': B_sky_BSFG < B_classical,
            'primary_equations': [
                'A_uv = g_uv + eta*T_s00*cos(pi*t_n)*delta_uv    [BSFG PAPER_554]',
                'R^r_0r0 = 6*eta*C_num/r^5',
                'E(r) = E_0*exp(-Gamma_BSFG*r)',
                'Gamma_BSFG = eta*|R_scalar_avg|/c^4',
                'B_n^BSFG = n*L*dr/(4*pi*c)*exp(-Gamma*r_n)*[SSq]^n',
                'B_sky^BSFG = sum_{n=1}^{26} B_n^BSFG  [VDS x BSFG double supp]',
            ],
            'available_equations': [
                'Combined bound: B_sky <= B_classical * Li_{26}([SSq]) * exp(-Gamma*r_H)',
                'BSFG horizon r_h=(eta*C_num)^{1/3}=0.233*R_sun [PAPER_561]',
                'r_cross=0.36AU Bohr-Sommerfeld crossover [PAPER_562]',
            ],
            'simulation_set': [
                'Vary eta in [1e-24,1e-20] -> Gamma_BSFG linear in eta',
                'Compare 3 UQFF methods: DPM / VDS / BSFG shell-by-shell',
                'Gap-fill: add n_star(z) Madau-Dickinson (PAPER_567)',
                'BSFG horizon blinking: E(r) pulsates with cos(pi*t_n)',
            ],
            'session': self.SESSION, 'papers': [self.PAPER],
        }


__all__ = ["""

if ANCHOR_BEFORE_ALL not in content:
    print('ANCHOR NOT FOUND in CP4 - searching for close match...')
    idx = content.find("            'session': self.SESSION, 'papers': [self.PAPER],\n        }\n\n\n__all__")
    if idx == -1:
        idx = content.find("\n__all__ = [")
        if idx == -1:
            print('ERROR: cannot find __all__ in CP4')
            sys.exit(1)
        # Insert before __all__
        INSERT_POINT = content.rfind('\n\n', 0, idx) + 2
        new_content = content[:INSERT_POINT] + NEW_CLASSES.rsplit('__all__ = [', 1)[0] + content[INSERT_POINT:]
        path.write_text(new_content, encoding='utf-8')
        print('Inserted (before __all__)')
    else:
        new_content = content[:idx] + NEW_CLASSES + content[idx + len("            'session': self.SESSION, 'papers': [self.PAPER],\n        }\n\n\n__all__ = ["):]
        path.write_text(new_content, encoding='utf-8')
        print('Inserted (partial anchor)')
else:
    new_content = content.replace(ANCHOR_BEFORE_ALL, NEW_CLASSES, 1)
    path.write_text(new_content, encoding='utf-8')
    print('Inserted via full anchor')
