import math

NEW_CLASS = """

# =============================================================================
# SESSION 193b CLASSES (v5.50) -- grok_share_96da8158-f7c5.txt FINAL PASS
# PAPER_827 -- W_stellar Orion/Eagle + P_term Hydrogen Atom
# FILE FULLY TAPPED - 100% EXTRACTION COMPLETE
# =============================================================================


class WstellarPtermOrionEagleHydrogenAtomUQFFCalculator:
    \\"\\"\\"PAPER_827 #411 -- W_stellar stellar wind pressure + P_term atomic pressure correction.
    W_stellar=Mdot_wind*v_wind/(4*pi*r^2*rho_cloud) (Orion/Eagle + P_rad net balance);
    P_term=(P_ext*a0^3*Z^2)/(E_1/n^2) multiplicative atomic pressure correction;
    Together complete 100% extraction of grok_share_96da8158-f7c5.txt.
    Documents 27 (Hydrogen Atom), 34 (Orion), 36 (Eagle). PAPER_827.
    \\"\\"\\"
    def compute(self, dataset: dict) -> dict:
        G = 6.6743e-11
        H0 = dataset.get('H0', 2.269e-18)
        hbar = 1.0546e-34
        c = 2.998e8
        lam = 1.1e-52
        t_Hubble = 4.352e17
        mode = dataset.get('mode', 'orion')  # 'orion', 'eagle', 'hydrogen'
        if mode in ('orion', 'eagle'):
            M = dataset.get('M', 2e33 if mode == 'orion' else 5e33)
            r = dataset.get('r', 3.09e15 if mode == 'orion' else 3.086e16)
            B = dataset.get('B', 0.0)
            B_crit = dataset.get('B_crit', 4.4e13)
            z = dataset.get('z', 4e-4 if mode == 'orion' else 2.0e-3)
            t = dataset.get('t', t_Hubble)
            rho_cloud = dataset.get('rho_cloud', 3.34e-24 if mode == 'orion' else 5e-24)
            Mdot_wind = dataset.get('Mdot_wind', 2.52e16 if mode == 'orion' else 6.3e16)
            v_wind = dataset.get('v_wind', 2e6 if mode == 'orion' else 1.5e6)
            L_star = dataset.get('L_star', 7.7e31 if mode == 'orion' else 1.5e36)
            kappa = dataset.get('kappa', 0.01)
            dx = dataset.get('delta_x', 1e-10)
            dp = dataset.get('delta_p', 1e-24)
            psi_sq = dataset.get('psi_sq', 1.0)
            M_vis = dataset.get('M_vis', M)
            M_DM = dataset.get('M_DM', 0.0)
            delta_rho_rho = dataset.get('delta_rho_rho', 0.01)
            # W_stellar
            W_stellar = Mdot_wind * v_wind / (4 * math.pi * r**2 * rho_cloud) if rho_cloud > 0 else 0.0
            # P_rad
            P_rad = L_star / (4 * math.pi * r**2 * c * rho_cloud * kappa) if rho_cloud > 0 else 0.0
            net_wind_rad = W_stellar - P_rad
            # Core gravity
            g_newton = G * M / r**2
            g_expansion = g_newton * (1 + H0 * t)
            g_superc = g_expansion * (1 - B / B_crit)
            g_cosmo = lam * c**2 / 3
            uncertainty = hbar / math.sqrt(dx * dp) if dx > 0 and dp > 0 else 0.0
            g_quantum = uncertainty * psi_sq * (2 * math.pi / t_Hubble)
            dm_term = (M_vis + M_DM) * (delta_rho_rho + 3 * G * M / r**3)
            g_total = g_superc + g_cosmo + g_quantum + dm_term + net_wind_rad
            pillar_state = 'wind-compressed' if net_wind_rad > 0 else 'photo-evaporating'
            return {
                'g_total': g_total,
                'W_stellar': W_stellar,
                'P_rad': P_rad,
                'net_wind_rad': net_wind_rad,
                'pillar_state': pillar_state,
                'mode': mode,
                'primary_equations': [
                    f"W_stellar=Mdot_wind*v_wind/(4*pi*r^2*rho_cloud)={W_stellar:.3e} m/s^2",
                    f"P_rad=L_star/(4*pi*r^2*c*rho_cloud*kappa)={P_rad:.3e} m/s^2",
                    f"net=W_stellar-P_rad={net_wind_rad:.3e} [{pillar_state}]",
                    f"g_total={g_total:.4e} m/s^2",
                ],
                'available_equations': [
                    "W_stellar = Mdot_wind*v_wind / (4*pi*r^2*rho_cloud)",
                    "P_rad = L_star / (4*pi*r^2*c*rho_cloud*kappa)",
                    "net > 0: pillar compression; net < 0: photoevaporation",
                ],
                'simulation_set': ["sweep r ionization-front", "vary Mdot_wind/L_star ratio", "mode orion/eagle"],
            }
        else:  # hydrogen atom
            m_p = 1.673e-27
            m_e = 9.109e-31
            a0 = 5.292e-11
            E1_J = 13.6 * 1.6e-19  # 13.6 eV in J
            n = dataset.get('n', 1)
            Z = dataset.get('Z', 1)
            r = dataset.get('r', a0 / Z)
            P_ext = dataset.get('P_ext', 2.4e-6)  # ambient radiation pressure
            F_tech = dataset.get('F_tech', 0.0)
            dx = dataset.get('delta_x', a0)
            dp = dataset.get('delta_p', hbar / a0)
            psi_sq = dataset.get('psi_sq', 1.0)
            E_n = E1_J / n**2
            t = dataset.get('t', t_Hubble)
            B = dataset.get('B', 0.0)
            B_crit = dataset.get('B_crit', 4.4e13)
            # P_term
            P_term = P_ext * a0**3 * Z**2 * n**2 / E1_J if E1_J > 0 else 0.0
            m_tot = m_p + m_e
            uncertainty = hbar / math.sqrt(dx * dp) if dx > 0 and dp > 0 else 0.0
            q_correction = uncertainty * psi_sq / E_n if E_n > 0 else 0.0
            g_atom = G * m_tot / r**2 * (1 + H0 * t) * (1 + P_term) * (1 + q_correction)
            g_cosmo = lam * c**2 / 3
            g_total = g_atom + g_cosmo + F_tech
            return {
                'g_total': g_total,
                'P_term': P_term,
                'F_tech': F_tech,
                'g_atom': g_atom,
                'mode': mode,
                'primary_equations': [
                    f"g_H=(G*(m_p+m_e))/r^2*(1+H0*t)*(1+P_term)*(1+q_corr)+...+F_tech",
                    f"P_term=P_ext*a0^3*Z^2*n^2/E_1={P_term:.3e} (dimensionless)",
                    f"g_total={g_total:.4e} m/s^2",
                ],
                'available_equations': [
                    "P_term = P_ext * a0^3 * Z^2 * n^2 / E_1",
                    "P_term significant when P_ext > ~1e10 Pa (laser/high-field lab)",
                    "F_tech = applied external field equivalent",
                ],
                'simulation_set': ["sweep P_ext 1e-6..1e12", "sweep n 1..7", "sweep Z 1..92"],
            }

    def simulate(self, sweep=None, sweep_param=None):
        modes = sweep or ['orion', 'eagle', 'hydrogen']
        results = []
        for m in modes:
            r = self.compute({'mode': m})
            r['sweep_val'] = m
            results.append(r)
        return results


# =============================================================================
# __all__ registry -- Session 193b additions
# =============================================================================

_SESSION_193b_CLASSES = [
    "WstellarPtermOrionEagleHydrogenAtomUQFFCalculator",
]
"""

with open('CondensedPhysics4.py', 'a', encoding='utf-8') as f:
    f.write(NEW_CLASS)

print('Done. Appended Session 193b class #411.')
