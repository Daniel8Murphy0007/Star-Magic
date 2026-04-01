"""
Session 173 — Generate _append_cp4_NNN.py files for PAPER_674-687, CP4 entries #258-271.
"""
import os, textwrap

ROOT = r"C:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic"
CP4  = os.path.join(ROOT, "CondensedPhysics2.py")

MODULES = [
    # (entry, PAPER, cp_class_name, cpp_module_name, description, primary_out, primary_in,
    #  [(param_name, type, default, desc), ...])
    (258, "PAPER_674", "UQFFComparedToLIGODataCalculator", "UQFFComparedToLIGOData",
     "UQFF vs General LIGO Dataset — chirp mass, h_UQFF(f), phase shift, frequency sweep",
     "h_uqff", "frequency_hz",
     [("m1_kg","float",7.16e31,"primary BH mass kg"),
      ("m2_kg","float",5.77e31,"secondary BH mass kg"),
      ("distance_m","float",1.27e25,"luminosity distance m"),
      ("frequency_hz","float",150.0,"GW frequency Hz")]),

    (259, "PAPER_675", "UQFFComparedToGW170817Calculator", "UQFFComparedToGW170817",
     "UQFF vs GW170817 NS-NS merger — GRB delay, tidal deformability, inspiral strain",
     "h_uqff", "frequency_hz",
     [("m1_kg","float",2.704e30,"NS1 mass kg (1.36 Msun)"),
      ("m2_kg","float",2.327e30,"NS2 mass kg (1.17 Msun)"),
      ("distance_m","float",1.234e24,"40 Mpc in m"),
      ("frequency_hz","float",1500.0,"GW frequency Hz")]),

    (260, "PAPER_676", "UQFFComparedToGW190425Calculator", "UQFFComparedToGW190425",
     "UQFF vs GW190425 heavy NS-NS — post-merger phase, ejecta limit, strain",
     "h_uqff", "frequency_hz",
     [("m1_kg","float",3.778e30,"NS1 mass kg (1.9 Msun)"),
      ("m2_kg","float",2.984e30,"NS2 mass kg (1.5 Msun)"),
      ("distance_m","float",4.9e24,"159 Mpc in m"),
      ("frequency_hz","float",2500.0,"GW frequency Hz")]),

    (261, "PAPER_677", "UQFFPredictionsForLISACalculator", "UQFFPredictionsForLISA",
     "UQFF predictions for LISA space GW observatory — h_LISA_UQFF, Omega_GW, EMRI rate",
     "h_uqff_lisa", "frequency_hz",
     [("chirp_mass_kg","float",1.989e38,"SMBH chirp mass (1e8 Msun)"),
      ("distance_m","float",1.0e27,"source distance m"),
      ("frequency_hz","float",1.0e-3,"GW frequency in LISA band Hz"),
      ("T_eff_K","float",2.73,"effective temperature K")]),

    (262, "PAPER_678", "LISAVsLIGOComparisonsCalculator", "LISAVsLIGOComparisons",
     "UQFF LISA vs LIGO suppression comparison — R_supp, crossover frequency, Sn_UQFF",
     "R_supp_ligo", "frequency_hz",
     [("mass_kg","float",3.978e31,"BH mass kg (20 Msun)"),
      ("frequency_hz","float",150.0,"GW frequency Hz"),
      ("T_eff_K","float",2.73,"effective temperature K")]),

    (263, "PAPER_679", "AetherSuperfluidDynamicsCalculator", "AetherSuperfluidDynamics",
     "UQFF aether superfluid — m_UA, healing length, sound speed, vortex circulation, g_eff",
     "g_eff", "radius_m",
     [("radius_m","float",1.0e10,"radial distance from BH m"),
      ("mass_kg","float",8.548e36,"BH mass kg (Sgr A*)"),
      ("n_UA","float",2.95e31,"aether number density m^-3")]),

    (264, "PAPER_680", "VortexQuantizationCalculator", "VortexQuantization",
     "UQFF vortex quantization — circulation kappa_v, core radius a_v, energy E_v, Magnus force",
     "vortex_energy", "winding_number",
     [("winding_number","int",1,"vortex winding number n"),
      ("R_outer_m","float",1.0e6,"outer boundary radius m"),
      ("v_s_m_s","float",1.0e3,"superfluid velocity m/s")]),

    (265, "PAPER_681", "GrossPitaevskiiVortexSimulationCalculator", "GrossPitaevskiiVortexSimulation",
     "UQFF GP vortex simulation — imaginary-time ground state, chemical potential, density profile",
     "chemical_potential_eV", "radius_m",
     [("r_min_m","float",1.0e9,"inner grid radius m"),
      ("r_max_m","float",1.0e11,"outer grid radius m"),
      ("mass_kg","float",8.548e36,"BH mass kg"),
      ("n_grid","int",100,"radial grid points"),
      ("n_steps","int",50,"imaginary-time steps")]),

    (266, "PAPER_682", "UQFFStabilityNumericallyForSgrACalculator", "UQFFStabilityNumericallyForSgrA",
     "UQFF numerical stability for Sgr A* — Lyapunov, omega_I, RK4 mass evolution, stability check",
     "lyapunov_exponent", "mass_kg",
     [("mass_kg","float",8.548e36,"BH mass kg (Sgr A*)"),
      ("t_end_s","float",1.0e60,"evolution end time s"),
      ("dt_s","float",1.0e55,"time step s")]),

    (267, "PAPER_683", "UQFFHawkingTemperatureModulationCalculator", "UQFFHawkingTemperatureModulation",
     "UQFF Hawking T modulation — T_UQFF, modulation factor, Wien peak, Planck spectrum shift",
     "T_uqff_K", "mass_kg",
     [("mass_kg","float",1.0e12,"BH mass kg"),
      ("radius_m","float",0.0,"evaluation radius m (0=r_s)"),
      ("time_s","float",1.0e8,"evaluation time s")]),

    (268, "PAPER_684", "UQFFPrimordialBHEvaporationCalculator", "UQFFPrimordialBHEvaporation",
     "UQFF primordial BH evaporation — dM/dt suppressed, tau_UQFF, M(t_form), RK4 evolution",
     "tau_uqff_s", "mass_kg",
     [("mass_kg","float",1.0e12,"PBH mass kg"),
      ("t_form_s","float",1.0e-23,"formation time s"),
      ("dt_s","float",1.0e50,"RK4 time step s")]),

    (269, "PAPER_685", "UQFFPBHDarkMatterImplicationsCalculator", "UQFFPBHDarkMatterImplications",
     "UQFF PBH dark matter — M_crit shift, f_PBH boost, viable mass window scan",
     "f_pbh_boost", "mass_kg",
     [("mass_kg","float",1.0e12,"PBH mass kg"),
      ("n_pbh_m3","float",1.0e6,"PBH number density m^-3")]),

    (270, "PAPER_686", "UQFFModulationForM87Calculator", "UQFFModulationForM87",
     "UQFF modulation for M87* — T_UQFF, shadow_UQFF, jet power, ring brightness, T_acc",
     "T_uqff_K", "mass_kg",
     [("mass_kg","float",1.293e40,"M87* mass kg (6.5e9 Msun)"),
      ("T_acc_K","float",1.0e7,"accretion disk temperature K")]),

    (271, "PAPER_687", "M87MassEvolutionSimulationCalculator", "M87MassEvolutionSimulation",
     "M87* mass evolution — Bondi-UQFF accretion, Hawking evaporation, jet mass loss, Hubble-time RK4",
     "dM_dt_total_kg_s", "mass_kg",
     [("mass_kg","float",1.293e40,"M87* mass kg"),
      ("rho_ism_kg_m3","float",1.67e-25,"ISM density kg/m^3"),
      ("T_ism_K","float",1.0e7,"ISM temperature K"),
      ("dt_s","float",1.0e13,"RK4 time step s")]),
]

def build_cp4_class(entry, paper, cp_cls, cpp_cls, description, primary_out, primary_in, params):
    """Generate a CondensedPhysics2.py calculator class source string."""
    param_lines = "\n".join(
        f"        ('{p[0]}', {repr(p[1])}, {repr(p[2])}, {repr(p[3])}),"
        for p in params
    )
    default_dict = ", ".join(
        f"'{p[0]}': {repr(p[2])}" for p in params
    )
    param_extract = "\n".join(
        f"        {p[0]} = params.get('{p[0]}', {repr(p[2])})" for p in params
    )
    compute_args = ", ".join(p[0] for p in params)

    # Physics implementation block
    phys_block = f"""\
        import math
        G      = 6.6743e-11
        C      = 2.998e8
        HBAR   = 1.0546e-34
        K_B    = 1.380649e-23
        PI     = math.pi
        M_SUN  = 1.989e30
        RHO_UA = 7.09e-36
        RHO_SCM= 7.09e-37
        F_TRZ  = 0.1
        KAPPA  = 0.0005
        SSQ    = 0.57
        MU_J   = 3.38e23
        GAMMA  = 5.0e-5 / 86400.0

        # Standard helpers
        def T_H(M):
            return HBAR*C**3 / (8.0*PI*G*M*K_B) if M>0 else 1e-300
        def L_H(M):
            return HBAR*C**6 / (15360.0*PI*G**2*M**2) if M>0 else 0
        def r_s(M):
            return 2.0*G*M / C**2 if M>0 else 0
        def tau_std(M):
            return 5120.0*PI*G**2*M**3 / (HBAR*C**4) if M>0 else 0
        def U_m(r, t):
            tn = t / 1.0e8
            return (MU_J/r)*(1.0 - math.exp(-GAMMA*t*math.cos(PI*tn))) if r>0 else 0

        # Common suppression factors
        def S_SCm(M_val, T_Hval, r_val):
            arg = RHO_SCM * r_val / (K_B * T_Hval) if T_Hval>0 else 0
            return math.exp(-min(arg, 700))
        def S_Um_f(M_val, r_val, t_val, f_val):
            Um = U_m(r_val, t_val)
            arg = Um * 2.0*PI*f_val / C**2
            return math.exp(-min(arg, 700))
"""

    # Per-class compute core
    paper_num = int(paper.split("_")[1])

    if paper_num == 674:
        compute_core = f"""\
        def chirp_mass(m1, m2):
            return (m1*m2)**0.6 / (m1+m2)**0.2
        def h_GR(Mc, d, f):
            return (4.0/d) * (G*Mc/C**2)**(5.0/3.0) * (PI*f)**(2.0/3.0) / C**2

        M  = m1_kg + m2_kg
        Mc = chirp_mass(m1_kg, m2_kg)
        r  = r_s(M); TH = T_H(M)
        S1 = S_SCm(M, TH, r)
        S2 = S_Um_f(M, r, 1e8, frequency_hz)
        h_gr   = h_GR(Mc, distance_m, frequency_hz)
        h_uqff = (1.0 - F_TRZ) * S1 * S2 * h_gr
        delta_phi = KAPPA * F_TRZ * 1e8
        result = {{'h_gr': h_gr, 'h_uqff': h_uqff, 'S_SCm': S1, 'S_Um': S2,
                   'chirp_mass_kg': Mc, 'delta_phi_rad': delta_phi}}
"""
    elif paper_num == 675:
        compute_core = f"""\
        def chirp_mass(m1, m2):
            return (m1*m2)**0.6 / (m1+m2)**0.2
        def h_GR(Mc, d, f):
            return (4.0/d) * (G*Mc/C**2)**(5.0/3.0) * (PI*f)**(2.0/3.0) / C**2
        M  = m1_kg + m2_kg; Mc = chirp_mass(m1_kg, m2_kg)
        r  = r_s(M); TH = T_H(M)
        S1 = S_SCm(M, TH, r)
        S2 = S_Um_f(M, r, 1e8, frequency_hz)
        h_gr   = h_GR(Mc, distance_m, frequency_hz)
        h_uqff = (1.0-F_TRZ)*S1*S2*h_gr
        grb_delay = 1.7*(1.0 + F_TRZ*(RHO_UA/RHO_SCM))
        kappa_tidal_ratio = 1.0 - F_TRZ*RHO_SCM/RHO_UA
        result = {{'h_gr': h_gr, 'h_uqff': h_uqff, 'grb_delay_UQFF_s': grb_delay,
                   'kappa_tidal_ratio': kappa_tidal_ratio}}
"""
    elif paper_num == 676:
        compute_core = f"""\
        def chirp_mass(m1, m2):
            return (m1*m2)**0.6 / (m1+m2)**0.2
        def h_GR(Mc, d, f):
            return (4.0/d) * (G*Mc/C**2)**(5.0/3.0) * (PI*f)**(2.0/3.0) / C**2
        M  = m1_kg + m2_kg; Mc = chirp_mass(m1_kg, m2_kg)
        r  = r_s(M); TH = T_H(M)
        S1 = S_SCm(M, TH, r); S2 = S_Um_f(M, r, 1e8, frequency_hz)
        h_gr   = h_GR(Mc, distance_m, frequency_hz)
        h_uqff = (1.0-F_TRZ)*S1*S2*h_gr
        ejecta_lim = 0.05*M * (RHO_SCM/RHO_UA)*(1.0-F_TRZ)
        result = {{'h_gr': h_gr, 'h_uqff': h_uqff, 'ejecta_limit_kg': ejecta_lim,
                   'post_merger_phase_rad': KAPPA*F_TRZ*1e-3}}
"""
    elif paper_num == 677:
        compute_core = f"""\
        def h_GR_SMBH(Mc, d, f):
            return (4.0/d) * (G*Mc/C**2)**(5.0/3.0) * (PI*f)**(2.0/3.0) / C**2
        L_LISA = 2.5e9
        S_UA   = max(0.0, 1.0 - RHO_UA*L_LISA/(K_B*T_eff_K))
        M_tot  = chirp_mass_kg * 4.0**0.2
        r  = r_s(M_tot); TH = T_H(M_tot)
        S1 = S_SCm(M_tot, TH, r)
        h_gr   = h_GR_SMBH(chirp_mass_kg, distance_m, frequency_hz)
        h_uqff_lisa = (1.0-F_TRZ)*S_UA*S1*h_gr
        omega_GW_ratio = (RHO_UA/9.47e-27)**F_TRZ
        EMRI_boost = 1.0 + F_TRZ*(RHO_UA/RHO_SCM)
        result = {{'h_uqff_lisa': h_uqff_lisa, 'S_UA_LISA': S_UA,
                   'omega_GW_UQFF_ratio': omega_GW_ratio, 'EMRI_rate_boost': EMRI_boost}}
        def chirp_mass(m1, m2): return (m1*m2)**0.6/(m1+m2)**0.2
"""
    elif paper_num == 678:
        compute_core = f"""\
        L_LISA = 2.5e9
        r  = r_s(mass_kg); TH = T_H(mass_kg)
        S1 = S_SCm(mass_kg, TH, r):
        S2 = S_Um_f(mass_kg, r, 1e8, frequency_hz)
        R_supp_ligo = (1.0-F_TRZ)*S1*S2
        S_UA = max(0.0, 1.0 - RHO_UA*L_LISA/(K_B*T_eff_K))
        R_supp_lisa = (1.0-F_TRZ)*S_UA*S1
        Um = U_m(r, 1e8)
        S_UA2 = max(1e-300, 1.0 - RHO_UA*L_LISA/(K_B*T_eff_K))
        if Um>0 and S_UA2>0:
            f_cross = -math.log(S_UA2)*C**2/(2.0*PI*Um)
        else:
            f_cross = 0.0
        result = {{'R_supp_ligo': R_supp_ligo, 'R_supp_lisa': R_supp_lisa,
                   'crossover_freq_Hz': f_cross}}
"""
    elif paper_num == 679:
        compute_core = f"""\
        H0_SI = 2.184e-18
        M_UA  = HBAR*H0_SI/C**2
        G_UA  = 1.0e-10
        N_UA  = RHO_UA/(HBAR*H0_SI/C**2)
        xi_UA = HBAR/math.sqrt(2.0*M_UA*G_UA*N_UA)
        c_UA  = math.sqrt(G_UA*N_UA/M_UA)
        kv    = 2.0*PI*HBAR/M_UA
        g_Newton = G*mass_kg/radius_m**2
        boost = 1.0+(c_UA**2/C**2)*F_TRZ*(RHO_UA/RHO_SCM)
        g_eff = g_Newton*boost
        result = {{'healing_length_m': xi_UA, 'sound_speed_m_s': c_UA,
                   'vortex_circulation': kv, 'g_eff': g_eff,
                   'g_Newton': g_Newton, 'm_UA_kg': M_UA}}
"""
    elif paper_num == 680:
        compute_core = f"""\
        H0_SI = 2.184e-18
        M_UA  = HBAR*H0_SI/C**2
        G_UA  = 1.0e-10
        N_UA  = RHO_UA/(HBAR*H0_SI/C**2)
        n = int(winding_number)
        kv    = n * 2.0*PI*HBAR/M_UA
        xi    = HBAR/math.sqrt(2.0*M_UA*G_UA*N_UA)
        a_v   = xi*math.exp(-n*PI)
        if a_v>0 and R_outer_m>a_v:
            Ev = RHO_UA*kv**2/(4.0*PI)*math.log(R_outer_m/a_v)*(RHO_UA/RHO_SCM)
        else:
            Ev = 0.0
        c_UA  = math.sqrt(G_UA*N_UA/M_UA)
        omega_v_dummy = c_UA/R_outer_m  # rough estimate
        omega_v_UQFF = omega_v_dummy*(1.0+F_TRZ*c_UA/C)
        F_Magnus = RHO_UA*kv*v_s_m_s*(1.0+F_TRZ*(RHO_UA/RHO_SCM))
        vortex_energy = Ev
        result = {{'vortex_energy': vortex_energy, 'circulation_kv': kv,
                   'core_radius_m': a_v, 'magnus_force_UQFF': F_Magnus,
                   'omega_v_UQFF': omega_v_UQFF}}
"""
    elif paper_num == 681:
        compute_core = f"""\
        H0_SI = 2.184e-18
        M_UA  = HBAR*H0_SI/C**2
        G_UA  = 1.0e-10
        N_UA  = RHO_UA/(HBAR*H0_SI/C**2)
        # Thomas-Fermi approx: mu = G_UA*n_UA
        mu_UA = G_UA*N_UA
        xi    = HBAR/math.sqrt(2.0*M_UA*G_UA*N_UA)
        c_UA  = math.sqrt(G_UA*N_UA/M_UA)
        # Gravitational potential at r_min_m
        V_grav = -G*mass_kg*M_UA/r_min_m
        E_total = mu_UA + V_grav
        chemical_potential_eV = E_total / 1.602e-19
        result = {{'chemical_potential_eV': chemical_potential_eV,
                   'healing_length_m': xi, 'mu_UA_J': mu_UA,
                   'TF_density_m3': N_UA, 'c_UA_m_s': c_UA}}
"""
    elif paper_num == 682:
        compute_core = f"""\
        M = mass_kg
        TH = T_H(M); r = r_s(M); tau = tau_std(M)
        Um = U_m(r, 1e8)
        omega_I = -1.0/max(tau,1e-300) * (1.0+F_TRZ*(RHO_UA/RHO_SCM)*Um/(K_B*TH))
        arg = Um/(K_B*TH) if TH>0 else 0
        lam = -(RHO_SCM/RHO_UA)*math.exp(-min(arg,700))/max(tau,1e-300)
        stable = lam < 0
        lyapunov_exponent = lam
        result = {{'lyapunov_exponent': lyapunov_exponent, 'omega_I': omega_I,
                   'stable': stable, 'T_H_K': TH, 'tau_std_s': tau}}
"""
    elif paper_num == 683:
        compute_core = f"""\
        M = mass_kg
        r = r_s(M) if radius_m<=0 else radius_m
        TH = T_H(M)
        Um = U_m(r, time_s)
        fac = (1.0+F_TRZ)*(1.0-RHO_SCM/RHO_UA)*(1.0+Um/(K_B*TH)) if TH>0 else 1.0
        T_uqff_K = TH * fac
        wien_GR   = HBAR*C/(2.82*K_B*TH) if TH>0 else 0
        wien_UQFF = HBAR*C/(2.82*K_B*T_uqff_K) if T_uqff_K>0 else 0
        result = {{'T_uqff_K': T_uqff_K, 'T_H_K': TH, 'mod_factor': fac,
                   'wien_GR_m': wien_GR, 'wien_UQFF_m': wien_UQFF}}
"""
    elif paper_num == 684:
        compute_core = f"""\
        M = mass_kg
        TH = T_H(M); r = r_s(M); tau = tau_std(M)
        Um = U_m(r, 1e8)
        dM_std = -HBAR*C**4/(15360.0*PI*G**2*M**2) if M>0 else 0
        arg = Um/(K_B*TH) if TH>0 else 0
        dM_UQFF = dM_std*(1.0-F_TRZ)*(RHO_SCM/RHO_UA)*math.exp(-min(arg,700))
        tau_UQFF_s = tau/(1.0-F_TRZ)*(RHO_UA/RHO_SCM)*math.exp(min(arg,700))
        # M at formation from t_form
        r_f = 2.998e8*t_form_s/2.0
        rho_rad = 3.0/(32.0*PI*G*t_form_s**2)
        M_form = (4.0/3.0)*PI*r_f**3*rho_rad
        result = {{'tau_uqff_s': tau_UQFF_s, 'tau_std_s': tau,
                   'dM_dt_UQFF': dM_UQFF, 'M_formation_kg': M_form,
                   'suppression_factor': (1.0-F_TRZ)*RHO_SCM/RHO_UA}}
"""
    elif paper_num == 685:
        compute_core = f"""\
        T_age = 4.34e17
        M_crit_std  = (HBAR*C**4*T_age/(5120.0*PI*G**2))**(1.0/3.0)
        tau_ratio   = 1.0/(1.0-F_TRZ)*(RHO_UA/RHO_SCM)
        M_crit_UQFF = M_crit_std / tau_ratio**(1.0/3.0)
        f_pbh_boost = tau_ratio**(2.0/3.0)
        rho_PBH     = mass_kg * n_pbh_m3 * f_pbh_boost
        viable_std  = mass_kg > M_crit_std
        viable_UQFF = mass_kg > M_crit_UQFF
        result = {{'f_pbh_boost': f_pbh_boost, 'M_crit_std_kg': M_crit_std,
                   'M_crit_UQFF_kg': M_crit_UQFF, 'rho_PBH_UQFF': rho_PBH,
                   'viable_std': viable_std, 'viable_UQFF': viable_UQFF}}
"""
    elif paper_num == 686:
        compute_core = f"""\
        M = mass_kg
        TH    = T_H(M)
        T_uqff_K = TH*(1.0+F_TRZ)*(1.0-RHO_SCM/RHO_UA)
        r_sh  = 3.0*math.sqrt(3.0)*G*M/C**2
        r_sh_UQFF = r_sh*math.sqrt(1.0+F_TRZ*(RHO_UA/RHO_SCM))
        P_JET_GR = 1.0e37
        P_jet_UQFF = P_JET_GR*(1.0+F_TRZ)*math.sqrt(RHO_UA/RHO_SCM)
        brightness_ratio = (RHO_UA/RHO_SCM)**(F_TRZ/4.0)
        r = r_s(M); Um = U_m(r, 1e8)
        T_acc_UQFF = T_acc_K*(1.0+Um/(K_B*TH)) if TH>0 else T_acc_K
        result = {{'T_uqff_K': T_uqff_K, 'T_H_K': TH,
                   'shadow_GR_m': r_sh, 'shadow_UQFF_m': r_sh_UQFF,
                   'jet_power_UQFF_W': P_jet_UQFF, 'brightness_ratio': brightness_ratio,
                   'T_acc_UQFF_K': T_acc_UQFF}}
"""
    elif paper_num == 687:
        compute_core = f"""\
        M = mass_kg
        # Bondi-UQFF accretion
        c_s = math.sqrt(5.0/3.0*K_B*T_ism_K/1.673e-27)
        Mdot_Bondi = 4.0*PI*G**2*M**2*rho_ism_kg_m3/c_s**3
        rho_eff = rho_ism_kg_m3 + RHO_UA - RHO_SCM
        Mdot_UQFF = Mdot_Bondi*(rho_eff/rho_ism_kg_m3)*(1.0+F_TRZ)
        # Evaporation
        TH = T_H(M); r = r_s(M); Um = U_m(r, 1e8)
        dM_std = -HBAR*C**4/(15360.0*PI*G**2*M**2) if M>0 else 0
        arg = Um/(K_B*TH) if TH>0 else 0
        dM_evap = dM_std*(1.0-F_TRZ)*(RHO_SCM/RHO_UA)*math.exp(-min(arg,700))
        # Jet loss
        rho_disk = rho_ism_kg_m3*1e6
        B_eq = math.sqrt(8.0*PI*rho_disk*C**2)
        if B_eq>1e8: B_eq=1e8
        Omega_H = C/(4.0*G*M/C**2)*C if M>0 else 0
        rs_m = r_s(M)
        P_BZ = 0.01*B_eq**2*rs_m**2*C*Omega_H**2/(4.0*PI*C)
        P_jet_UQFF = P_BZ*(1.0+F_TRZ)*math.sqrt(RHO_UA/RHO_SCM)
        dM_jet = -P_jet_UQFF/C**2
        dM_dt_total_kg_s = Mdot_UQFF + dM_evap + dM_jet
        result = {{'dM_dt_total_kg_s': dM_dt_total_kg_s, 'Mdot_Bondi_UQFF': Mdot_UQFF,
                   'dM_evap': dM_evap, 'dM_jet': dM_jet, 'P_jet_W': P_jet_UQFF}}
"""
    else:
        compute_core = f"        result = {{}}\n"

    # Fix 678 syntax error ('S1 = S_SCm(...):' has colon)
    if paper_num == 678:
        compute_core = compute_core.replace(
            "        S1 = S_SCm(mass_kg, TH, r):",
            "        S1 = S_SCm(mass_kg, TH, r)"
        )

    class_body = f"""\
class {cp_cls}:
    \"\"\"
    CP4 #{entry} — {paper} | {description}
    CPP_MODULE: {cpp_cls}
    \"\"\"
    ENTRY   = {entry}
    PAPER   = "{paper}"
    CPP_MODULE = "{cpp_cls}"
    UQFF_CONSTANTS = {{
        "RHO_UA":  7.09e-36,
        "RHO_SCM": 7.09e-37,
        "F_TRZ":   0.1,
        "KAPPA":   0.0005,
        "SSQ":     0.57,
        "MU_J":    3.38e23,
        "GAMMA":   5.0e-5 / 86400.0,
    }}
    PARAMETERS = [
{param_lines}
    ]
    PRIMARY_OUTPUT = "{primary_out}"
    PRIMARY_INPUT  = "{primary_in}"

    def compute(self, params: dict = None) -> dict:
        if params is None:
            params = {{}}
        defaults = {{{default_dict}}}
        for k, v in defaults.items():
            params.setdefault(k, v)
{param_extract}
{compute_core}        self._last = result
        return result

    def simulate(self, sweep: list = None, sweep_param: str = None) -> list:
        if sweep is None:
            sweep = [{{}}]
        results = []
        for p in sweep:
            results.append(self.compute(dict(p)))
        return results

    def add_mod(self, fn) -> None:
        if not hasattr(self, "_mods"):
            self._mods = []
        self._mods.append(fn)

    def update_from_file(self, filepath: str) -> None:
        import os
        if not os.path.isfile(filepath):
            return
        with open(filepath, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if "=" in line:
                    k, v = line.split("=", 1)
                    try:
                        object.__setattr__(self, k.strip(), float(v.strip()))
                    except ValueError:
                        pass
"""
    return class_body


def write_append_script(entry, paper, cp_cls, class_src):
    script_path = os.path.join(ROOT, f"_append_cp4_{entry}.py")
    divider = (
        f"\n\n# {'='*72}\n"
        f"# CP4 #{entry} \u2014 {cp_cls}\n"
        f"# {paper} | Session 173 | appended by _append_cp4_{entry}.py\n"
        f"# {'='*72}\n"
    )
    class_src_escaped = class_src.replace("\\", "\\\\").replace('"""', '\\"\\"\\"')
    script = f"""#!/usr/bin/env python3
\"\"\"Append CP4 #{entry} ({cp_cls}) to CondensedPhysics2.py\"\"\"
import os, sys

CP4  = r"{CP4}"
DIVIDER = {repr(divider)}
CLASS_SRC = {repr(class_src)}

def main():
    if not os.path.isfile(CP4):
        print(f"ERROR: Cannot find {{CP4}}", file=sys.stderr); sys.exit(1)
    # Check if already appended
    with open(CP4, "rb") as f:
        existing = f.read().decode("utf-8", errors="replace")
    if "class {cp_cls}" in existing:
        print(f"SKIP: {cp_cls} already present in CondensedPhysics2.py")
        return
    block = DIVIDER + CLASS_SRC
    with open(CP4, "ab") as f:
        f.write(block.encode("utf-8"))
    print(f"APPENDED: CP4 #{entry} {cp_cls} -> CondensedPhysics2.py")
    # Verify
    with open(CP4, "rb") as f:
        tail = f.read()[-2000:].decode("utf-8", errors="replace")
    if "class {cp_cls}" in tail or "{cp_cls}" in existing:
        print(f"VERIFIED: {cp_cls} present")
    else:
        print(f"WARNING: Could not verify {cp_cls} — check file")

if __name__ == "__main__":
    main()
"""
    with open(script_path, "w", encoding="utf-8") as f:
        f.write(script)
    print(f"  CREATED: {os.path.basename(script_path)}")


if __name__ == "__main__":
    os.chdir(ROOT)
    print("Generating CP4 append scripts for PAPER_674-687...")
    for (entry, paper, cp_cls, cpp_cls, desc, p_out, p_in, params) in MODULES:
        class_src = build_cp4_class(entry, paper, cp_cls, cpp_cls, desc, p_out, p_in, params)
        write_append_script(entry, paper, cp_cls, class_src)
    print(f"Done: 14 _append_cp4_NNN.py files created.")
