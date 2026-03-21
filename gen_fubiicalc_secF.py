"""
gen_fubiicalc_secF.py
Returns C++ buildSectionF() — 68 unique Um (Universal Magnetism) variant types.
Source: BB_C_Equations_04Sept2025.pdf (extracted from grok_share_c020496d9e.txt ~lines 693-2100).

General form for all Um variants:
  Um,TYPE = Sum_j [mu_j(t,rho_vac[SCm])/r_j * (1-exp(-gamma*t)*cos(pi*t_n)) * phi^j]
          * P_SCm * E_react * (1+1e13*f_Heaviside) * (1+f_quasi) * exp(-[SSq])
          * TYPE_SPECIFIC_FACTOR
"""

# (id, name, expression, context, validation_note, numerical_result, units, system)
_F_DATA = [
    (693,  "Um_general",
     "Um = Sum_j[mu_j(t,rho_vac[SCm])/r_j*(1-exp(-gamma*t)*cos(pi*t_n))*phi^j]"
     "*P_SCm*E_react*(1+1e13*f_Heaviside)*(1+f_quasi)*exp(-[SSq])",
     "Universal Magnetism general form - all-scale vacuum magnetic moment summation",
     "gamma~5e-5 day^-1; phi~1.02; Um~3.78e-6 J/m^3; base of all 68 variants",
     3.78e-6, "J/m^3", "Universal magnetism"),

    (724,  "Um_cosmo",
     "Um_cosmo = Um_general * (k*c^2/a^2) / H^2",
     "Cosmological Um - scaled by curvature/Hubble ratio; FRW correction",
     "k=curvature; a=scale factor; H=Hubble parameter; applies to cosmic epochs",
     0.0, "J/m^3", "Cosmological Um"),

    (746,  "Um_reh",
     "Um_reh = Um_general * (30*V_end/(pi^2*g_star))^(1/4) * exp(-3*N_reh/4)",
     "Reheating Um - inflaton decay epoch magnetic moment; reheat temperature factor",
     "V_end=inflaton potential at end; g_star~106; N_reh e-folds of reheating",
     0.0, "J/m^3", "Reheating Um"),

    (748,  "Um_MHD",
     "Um_MHD = Um_general * 4*pi*rho*v_turb * M_A",
     "MHD turbulent Um - magnetic energy from turbulent velocity * Alfven Mach M_A",
     "M_A=v_turb/v_A; v_A=Alfven speed; sub-Alfvenic turbulence M_A<1",
     0.0, "J/m^3", "MHD turbulent Um"),

    (770,  "Um_BBN",
     "Um_BBN = Um_general * (3/(32*pi*G*rho_rad)); t_BBN ~ 180 s",
     "BBN epoch Um - cosmic time at nucleosynthesis; Friedmann time-density",
     "t_BBN~180 s at T~0.1 MeV; primordial B-field constraint from BBN",
     0.0, "J/m^3", "BBN Um"),

    (802,  "Um_qnm",
     "Um_qnm = Um_general * a*c^3/(2*G*M); l=2, m=2",
     "Quasi-normal mode Um - BH ringdown; Kerr spin a; dominant harmonic",
     "omega_QNM = c^3/(2piGM)*(0.3737+0.088*a); tau_QNM ~ M",
     0.0, "J/m^3", "QNM Um"),

    (804,  "Um_bz",
     "Um_bz = Um_general * B^2*Omega_H^2*R_H^4; ~ a^2",
     "Blandford-Znajek Um - jet power magnetic factor; horizon angular velocity",
     "Omega_H=a*c^3/(2*G*M*r_H); P_BZ~a^2 for small spin",
     0.0, "J/m^3", "BZ Um"),

    (830,  "Um_photo",
     "Um_photo = Um_general * eps_photo; eps_photo = 0.1 to 0.3 * L_X/F_X",
     "Photoevaporation Um - X-ray efficiency fraction; EUV/FUV heating",
     "xi=L_X/(n*r^2); eps~0.1-0.3; planetary atmosphere evaporation",
     0.0, "J/m^3", "Photoevaporation Um"),

    (832,  "Um_mig",
     "Um_mig = Um_general * Sigma*(H/r)^(-3)",
     "Migration Um - disk surface density / aspect ratio; Type I torque scaling",
     "Sigma=surface density; H/r~0.05; migration rate ~ Sigma*(H/r)^(-3)",
     0.0, "J/m^3", "Migration Um"),

    (854,  "Um_termv",
     "Um_termv = Um_general * (L/L_Edd); ~ 1",
     "Terminal velocity Um - Eddington luminosity ratio; stellar wind driving",
     "L/L_Edd~1 for OB stars; v_inf~3*v_esc; mass-loss rate critical",
     0.0, "J/m^3", "Terminal velocity Um"),

    (878,  "Um_ps",
     "Um_ps = Um_general * sigma_M_rms",
     "Press-Schechter Um - rms density fluctuation on scale M; mass variance",
     "sigma_M = rms of delta smoothed on sphere of mass M; sigma~1 at M*",
     0.0, "J/m^3", "Press-Schechter Um"),

    (880,  "Um_sfe",
     "Um_sfe = Um_general * mdot_halo; mdot ~ M^(1.1)*(1+z)^(2.25)",
     "Star formation efficiency Um - halo accretion rate scaling",
     "Dekel+2013; cold stream mass flux; SFR~mdot_halo*f_b*eps_sf",
     0.0, "J/m^3", "SFE Um"),

    (904,  "Um_bhent",
     "Um_bhent = Um_general; dM = T_H*dS",
     "BH entropy Um - thermodynamic 1st law; Hawking temperature and entropy",
     "dS_BH=dA/(4*l_pl^2); T_H*dS=dM*c^2; 1st law of BH mechanics",
     0.0, "J/m^3", "BH entropy Um"),

    (906,  "Um_evapl",
     "Um_evapl = Um_general; dM/dt = -P_evap/c^2; P~hbar*c^6/(G^2*M^2)",
     "Evaporation lifetime Um - BH Hawking radiation power rate",
     "P_Hawking~hbar*c^6/(15360*pi*G^2*M^2); tau~M^3; Planck mass remnant?",
     0.0, "J/m^3", "BH evaporation Um"),

    (924,  "Um_ang",
     "Um_ang = Um_general * G*M/r^3 * r_A",
     "Angular momentum Um - specific angular momentum at Alfven radius r_A",
     "r_A=Alfven radius in stellar wind; angular momentum shed per unit mass",
     0.0, "J/m^3", "Angular momentum Um"),

    (926,  "Um_jetvel",
     "Um_jetvel = Um_general * Omega*r_A*(r_A/r0)^(1/2)",
     "Jet velocity Um - magneto-centrifugal launch; Alfven lever arm",
     "Blandford-Payne 1982; r_A/r0=lever arm ratio; efficiency eta_MHD",
     0.0, "J/m^3", "Jet velocity Um"),

    (944,  "Um_chirp",
     "Um_chirp = Um_general * (32/5)*(G*M_c^(5/3)/c^5)*(pi*f)^(10/3)",
     "Chirp mass GW Um - inspiral power spectral density",
     "M_c=chirp mass; GW luminosity P_GW~M_c^(10/3)*f^(10/3)",
     0.0, "J/m^3", "Chirp mass Um"),

    (992,  "Um_pec",
     "Um_pec = Um_general * (-(f*H/(4*pi))*grad^-1*(delta))",
     "Peculiar velocity Um - linear density-velocity relation",
     "f=d ln D/d ln a~Omega_m^0.55; redshift space distortion",
     0.0, "J/m^3", "Peculiar velocity Um"),

    (994,  "Um_ion",
     "Um_ion = Um_general * Int(alpha_B*n_e^2*C*dt)",
     "Ionisation Um - recombination rate integral; clumping factor C",
     "alpha_B=case B recombination; C=<n^2>/<n>^2 ~3; reionization budget",
     0.0, "J/m^3", "Reionisation Um"),

    (1012, "Um_deb",
     "Um_deb = Um_general * rho_rad; rho_rad=pi^2*k_B^4*T^4/(30*hbar^3*c^5)*g_star",
     "Debye/radiation epoch Um - radiation energy density at temperature T",
     "g_star=106.75 SM; scales as T^4; dominates for T>Teq~3000 K",
     0.0, "J/m^3", "Radiation epoch Um"),

    (1014, "Um_fried1",
     "Um_fried1 = Um_general * sqrt(Omega_m*(1+z)^3+Omega_k*(1+z)^2+Omega_Lambda)",
     "Friedmann Um - H(z) normalised Hubble parameter; all energy components",
     "H(z)=H0*E(z); E^2=Omega_m*(1+z)^3+Omega_k*(1+z)^2+Omega_Lambda",
     0.0, "J/m^3", "Friedmann Um"),

    (1034, "Um_damp",
     "Um_damp = Um_general * Q_QNM; Q~10 for l=2, m=2",
     "QNM damping Um - quality factor Q of ringdown oscillation",
     "Q=pi*f*tau; l=2,m=2 most excited; Q~2-12 for a=0-0.99",
     0.0, "J/m^3", "QNM damping Um"),

    (1036, "Um_arnett",
     "Um_arnett = Um_general * t_d^2; t_d^2 = 3*kappa*M_ej/(4*pi*c*v_ej^2)",
     "Arnett diffusion time Um - SN Ia light curve peak timescale",
     "Arnett 1982; M_max ~ mdot_Ni; t_rise~20 days for typical SN Ia",
     0.0, "J/m^3", "Arnett SN diffusion Um"),

    (1058, "Um_dyn",
     "Um_dyn = Um_general * eta_resistivity",
     "Dynamo resistivity Um - Ohmic diffusivity eta; reconnection coupling",
     "eta=c^2/(4*pi*sigma); resistive MHD; fast reconnection for large Rm",
     0.0, "J/m^3", "Dynamo resistivity Um"),

    (1080, "Um_metal",
     "Um_metal = Um_general; d(Z*M_gas)/dt = Y*SFR - Z*mdot_out",
     "Metallicity Um - metal production minus outflow; yield Y; loading",
     "Effective yield Y_eff=Y/(1+eta_load); mass-metallicity scatter",
     0.0, "J/m^3", "Metallicity Um"),

    (1082, "Um_cool",
     "Um_cool = Um_general * n_e*n_i*Lambda_cool(T)",
     "Cooling Um - bremsstrahlung + line cooling; Lambda(T) cooling function",
     "Lambda~10^-23 erg cm^3/s at T~10^7 K; tcool~3kT/(n*Lambda)",
     0.0, "J/m^3", "Cooling function Um"),

    (1100, "Um_rev",
     "Um_rev = Um_general * eta_resistivity",
     "Reversal Um - same resistivity coupling as Um_dyn but in reversal context",
     "Magnetic polarity reversal; alpha-Omega dynamo; l_rev scale",
     0.0, "J/m^3", "Magnetic reversal Um"),

    (1102, "Um_wde",
     "Um_wde = Um_general * rho_DE; rho_DE ~ a^(-3*(1+w))",
     "w-dark energy Um - quintessence density evolution with w",
     "w=-1: constant; w>-1: decreasing; tracker solutions phi^2/2~V",
     0.0, "J/m^3", "w-dark energy Um"),

    (1252, "Um_vir",
     "Um_vir = Um_general; sigma_v^2 = G*M/(3*r) for gas; 2*KE+PE=0",
     "Virial Um - equilibrium cluster velocity dispersion",
     "T_vir=m*G*M/(3*k_B*r); X-ray temperature diagnostic",
     0.0, "J/m^3", "Virial Um"),

    (1272, "Um_ms",
     "Um_ms = Um_general * L_MS; L_MS ~ mu^4*M^3 [Eddington main seq.]",
     "Main sequence luminosity Um - mass-luminosity scaling",
     "mu=mean molecular weight; L~M^4 for stars; opacity-dominated",
     0.0, "J/m^3", "Main sequence Um"),

    (1274, "Um_ml",
     "Um_ml = Um_general; L ~ M^3.5 from HR diagram empirical fit",
     "Mass-luminosity Um - observational fit from HR diagram",
     "Empirical; L/Lsun~(M/Msun)^3.5; valid for 0.1-50 Msun",
     0.0, "J/m^3", "Mass-luminosity Um"),

    (753,  "Um_star",
     "Um_star = Um_general * (sigma_v^3/(G^2*rho*m*lnLambda))",
     "Star cluster Um - relaxation timescale coupling; stellar encounters",
     "t_re=0.34*sigma_v^3/(G^2*rho*m*lnLambda); globular cluster dynamics",
     0.0, "J/m^3", "Star cluster Um"),

    (881,  "Um_fb",
     "Um_fb = Um_general * eta_fb*mdot_star*c^2",
     "Feedback Um - efficiency-weighted energy injection rate",
     "eta_fb~0.001 stellar winds; ~0.1 AGN; quenches SF in massive halos",
     0.0, "J/m^3", "Feedback Um"),

    (899,  "Um_grow",
     "Um_grow = Um_general * delta; delta_ddot+2*H*delta_dot=(3/2)*Omega_m*H^2*delta",
     "Growth Um - density perturbation amplitude coupling",
     "D(a)~a matter dom; growth rate f=d ln D/d ln a~Omega_m^0.55",
     0.0, "J/m^3", "Growth Um"),

    (901,  "Um_haw",
     "Um_haw = Um_general * T_H; T_H=hbar*c^3/(8*pi*G*M*k_B)",
     "Hawking temperature Um - thermal radiation coupling from BH horizon",
     "T_H~6e-8 K (Msun); BH magnetic moment analogue in thermal field",
     0.0, "J/m^3", "Hawking Um"),

    (907,  "Um_lqcf",
     "Um_lqcf = Um_general * H_LQC; H^2=(8*pi*G*rho/3)*(1-rho/rho_crit)",
     "LQC Friedmann Um - loop quantum bounce Hubble rate coupling",
     "Polymerisation correction; maximum rho=rho_crit=0.41*rho_Pl",
     0.0, "J/m^3", "LQC Friedmann Um"),

    (911,  "Um_bounc",
     "Um_bounc = Um_general * t_bounce; t_bounce=sqrt(3/(16*pi*G*rho_crit))",
     "Bounce time Um - LQC quantum bounce timescale",
     "t_bounce~Planck time; connects pre/post-bounce epochs",
     0.0, "J/m^3", "LQC bounce Um"),

    (913,  "Um_ent",
     "Um_ent = Um_general * S_ent*exp(-t/tau_dec); S=-Tr(rho_A*ln rho_A)",
     "Entanglement entropy Um - decoherence-weighted von Neumann entropy",
     "tau_dec=decoherence time; AdS/CFT Ryu-Takayanagi formula",
     0.0, "J/m^3", "Entanglement Um"),

    (919,  "Um_jshock",
     "Um_jshock = Um_general * (rho1*v1^2+P1)",
     "J-shock Um - total momentum flux across discontinuous shock",
     "Rankine-Hugoniot; adiabatic index gamma=5/3 monatomic gas",
     0.0, "J/m^3", "J-shock Um"),

    (921,  "Um_cshock",
     "Um_cshock = Um_general * rho_n*v_n*nu_ni",
     "C-shock Um - ion-neutral friction rate; continuous shock structure",
     "nu_ni~1.5e-9 cm^3/s * n_i; L_d=v_di/nu_ni",
     0.0, "J/m^3", "C-shock Um"),

    (923,  "Um_halo",
     "Um_halo = Um_general * (dN_halo/dM_dz)",
     "Halo formation Um - merger tree rate coupling",
     "Lacey-Cole 1993; conditional probability; EPS framework",
     0.0, "J/m^3", "Halo Um"),

    (925,  "Um_disk",
     "Um_disk = Um_general * eps*M_gas/t_dyn",
     "Disk SF Um - gravitational instability star formation rate",
     "Q<1 for instability; Kennicutt-Schmidt; Sigma_SFR~Sigma_gas^1.4",
     0.0, "J/m^3", "Disk SF Um"),

    (927,  "Um_burst",
     "Um_burst = Um_general * mdot_gas_inflow*eps_burst",
     "Starburst Um - merger-driven gas inflow rate",
     "eps_burst~0.5-1; triggered by tidal compression; Arp 220 analogue",
     0.0, "J/m^3", "Starburst Um"),

    (929,  "Um_sedov",
     "Um_sedov = Um_general * (E*t^2/rho)^(1/5)",
     "Sedov-Taylor Um - blast wave radius coupling to magnetism",
     "R_sedov~t^(2/5); transitions to radiative at R~30 pc",
     0.0, "J/m^3", "Sedov-Taylor Um"),

    (931,  "Um_dsa",
     "Um_dsa = Um_general * u_s^2/c^2",
     "DSA Um - shock velocity squared coupling; particle acceleration",
     "N(E)~E^(-2) for strong shock; maximum E_max=Z*e*B*R_Hillas",
     0.0, "J/m^3", "DSA Um"),

    (941,  "Um_tov",
     "Um_tov = Um_general * (dP/dr_TOV)",
     "TOV Um - neutron star pressure gradient coupling",
     "Stiff EOS: R~13 km; soft EOS: R~10 km; GW170817 R>10.7 km",
     0.0, "J/m^3", "TOV Um"),

    (995,  "Um_acc",
     "Um_acc = Um_general * (L_Edd/(eps*c^2))*exp(t/t_Sal)",
     "Accretion Um - Eddington rate exponential growth; Salpeter time",
     "t_Sal=45 Myr; e-folding to SMBH from seed; quasar duty cycle",
     0.0, "J/m^3", "Eddington accretion Um"),

    (981,  "Um_IGM",
     "Um_IGM = Um_general * k_B*T_IGM/m_p",
     "IGM thermal Um - post-reionisation warm-hot gas coupling",
     "T_IGM~1e4-1e7 K; WHIM 30-40% of baryons; Ly-alpha forest",
     0.0, "J/m^3", "IGM Um"),

    (983,  "Um_gal",
     "Um_gal = Um_general * M_halo*H(z)^2",
     "Galaxy formation Um - halo mass times Hubble growth factor",
     "SFR~M_halo*f_b*eps_sf*H(z); cosmic SFR peaks z~2",
     0.0, "J/m^3", "Galaxy formation Um"),

    (985,  "Um_quant",
     "Um_quant = Um_general * H^2/(8*pi^2*eps*M_pl^2)",
     "Quantum inflation Um - inflationary power spectrum amplitude",
     "A_s~2.2e-9; eps~0.004 Starobinsky; pivot k_0=0.05 Mpc^-1",
     0.0, "J/m^3", "Quantum inflation Um"),

    (989,  "Um_w",
     "Um_w = Um_general * rho_DE*(1+w)",
     "DE equation-of-state Um - quintessence pressure coupling",
     "w+1 tracks field kinetic energy; w~-0.95 DESI 2024",
     0.0, "J/m^3", "DE EOS Um"),

    (991,  "Um_ng",
     "Um_ng = Um_general * f_NL; f_NL<10 (Planck 2018)",
     "Non-Gaussianity Um - primordial bispectrum coupling",
     "Local f_NL; equilateral f_NL; orthogonal f_NL; multi-field inflation",
     0.0, "J/m^3", "Non-Gaussianity Um"),

    (965,  "Um_nfwrot",
     "Um_nfwrot = Um_general * v_NFW^2; v^2=4*pi*G*rho_s*r_s^2*(ln(1+x)-x/(1+x))/r",
     "NFW rotation curve Um - halo profile velocity coupling",
     "rho_s*r_s^3=const; concentration c=r_vir/r_s; c~10-15 Milky Way",
     0.0, "J/m^3", "NFW Um"),

    (967,  "Um_sidm",
     "Um_sidm = Um_general * (sigma_m/m)",
     "SIDM cross-section Um - self-interaction velocity-dependent coupling",
     "sigma/m~1 cm^2/g for core formation; v^(-4) for light mediator",
     0.0, "J/m^3", "SIDM Um"),

    (772,  "Um_void",
     "Um_void = Um_general * (-3/5)*(Omega_m*a+Omega_Lambda)^(-3/2)*delta_v0",
     "Void density Um - linear underdense expansion",
     "Voids evacuate at rate ~ a^(-1); supervoids ~100 Mpc empty",
     0.0, "J/m^3", "Void Um"),

    (774,  "Um_reion",
     "Um_reion = Um_general * ndot_gamma*eps_esc",
     "Reionization Um - ionising photon emission rate coupling",
     "eps_esc~0.1-0.2; z_reion~6-10; completed by z~6 (JWST 2024)",
     0.0, "J/m^3", "Reionization Um"),

    (776,  "Um_ISM",
     "Um_ISM = Um_general * pi*c_s^2/(G*rho)",
     "ISM Jeans Um - thermal support scale coupling to magnetism",
     "lambda_J=sqrt(pi*c_s^2/(G*rho)); B-field adds magnetic Jeans scale",
     0.0, "J/m^3", "ISM Um"),

    (782,  "Um_merg2",
     "Um_merg2 = Um_general * (5/(256*pi^(8/3)))*c^5*M_c^(-5/3)*f^(-11/3)",
     "Merger coalescence Um - GW inspiral frequency-mass coupling",
     "Time to coalescence from f_i; dE/dt=-P_GW; orbit shrinks",
     0.0, "J/m^3", "GW coalescence Um"),

    (739,  "Um_CR",
     "Um_CR = Um_general * E_max * Z",
     "Cosmic ray Um - maximum energy times charge; Hillas criterion",
     "E_max=Z*e*B*R; Auger UHECRs E>1e19 eV; ankle E~3e18 eV",
     0.0, "J/m^3", "Cosmic ray Um"),

    (821,  "Um_agn",
     "Um_agn = Um_general * L_AGN/c",
     "AGN radiation Um - momentum flux from quasar luminosity",
     "L_AGN up to 10^48 erg/s; momentum coupling to host ISM",
     0.0, "J/m^3", "AGN Um"),

    (1965, "Um_jeans",
     "Um_jeans = Um_general * sqrt(pi*c_s^2/(G*rho))",
     "Jeans scale Um - instability wavelength magnetic coupling",
     "M_J=(4*pi/3)*(lambda_J/2)^3*rho; cloud collapse criterion",
     0.0, "J/m^3", "Jeans Um"),

    (1985, "Um_relax",
     "Um_relax = Um_general * t_relax; t_relax=0.34*sigma_v^3/(G^2*rho*m*lnLambda)",
     "Relaxation Um - two-body relaxation timescale coupling",
     "Phase mixing; ergodic distribution; violent relaxation Lynden-Bell",
     0.0, "J/m^3", "Relaxation Um"),

    (1945, "Um_voidden",
     "Um_voidden = Um_general * delta_v(a); delta_v~a^(-1)",
     "Void density contrast Um - underdense region evolution",
     "Linear density contrast in void; compensated void profile",
     0.0, "J/m^3", "Void density Um"),

    (2005, "Um_conv",
     "Um_conv = Um_general * v_conv; v_conv~(alpha^2*g*deltaT*H_p/(4*T))^(1/3)",
     "Convective Um - mixing-length velocity magnetic coupling",
     "alpha=MLT parameter; H_p=pressure scale height; stellar convection zones",
     0.0, "J/m^3", "Convective Um"),

    (1805, "Um_eta",
     "Um_eta = Um_general * eta_b; eta_b=n_b/n_gamma=6e-10",
     "Baryon-photon ratio Um - BBN era baryon asymmetry coupling",
     "Measured via D/H primordial; consistent with CMB Omega_b",
     6.0e-10, "J/m^3", "Baryon-photon Um"),

    (1846, "Um_GW",
     "Um_GW = Um_general * rho_GW; rho_GW=(c^2/(32*pi*G))*<hdot_ij^2>",
     "Gravitational wave Um - GW energy density coupling to magnetism",
     "PTA nHz background; LIGO HF; UQFF connects GW to Um",
     0.0, "J/m^3", "Gravitational wave Um"),

    (1847, "Um_inf",
     "Um_inf = Um_general * V(phi); V~(3*H^2*M_pl^2/8*pi)*exp(N)",
     "Inflation potential Um - slow-roll potential energy coupling",
     "Chaotic/Starobinsky/natural inflation; r<0.06 Planck",
     0.0, "J/m^3", "Inflation potential Um"),

    (1848, "Um_glitch",
     "Um_glitch = Um_general * Delta_Omega/Omega; Delta_Omega~1e-6 to 1e-4",
     "Pulsar glitch Um - fractional spin-up coupling to magnetism",
     "Vortex unpinning; crustal superfluid; Vela glitch Delta_nu/nu~1e-6",
     0.0, "J/m^3", "Pulsar glitch Um"),

    (1849, "Um_after",
     "Um_after = Um_general * F_nu_afterglow; F_nu ~ nu^(-(p-1)/2)*t^(-alpha_t)",
     "GRB afterglow Um - synchrotron flux spectrum coupling",
     "p~2.2; alpha_t=3(p-1)/4; jet break steepens decay",
     0.0, "J/m^3", "GRB afterglow Um"),

    (1850, "Um_cmb_pol",
     "Um_cmb_pol = Um_general * C_l_BB; CMB B-mode polarisation power",
     "CMB B-mode Um - primordial GW through lensing B-mode coupling",
     "r<0.06 Planck; BICEP/Keck; dusty CMB foreground subtraction challenge",
     0.0, "J/m^3", "CMB B-mode Um"),

    (1851, "Um_duty",
     "Um_duty = Um_general * (1-exp(-t/tau_cool))*(1+mdot_acc/mdot_Edd)^(-1)",
     "AGN duty cycle Um - time-averaged accretion activity coupling",
     "Duty_cycle~0.01-0.1; merger-triggered AGN episodes; calorimetry",
     0.0, "J/m^3", "AGN duty Um"),
]


def _esc(s: str) -> str:
    """Escape backslashes and double-quotes for C++ string literals."""
    return s.replace("\\", "\\\\").replace('"', '\\"')


def get_section_F() -> str:
    lines = [
        "\n// ======================================================\n",
        "// SECTION F:  68 unique Um (Universal Magnetism) variant types\n",
        "// Source: BB_C_Equations_04Sept2025.pdf\n",
        "//         (grok_share_c020496d9e.txt ~lines 693-2100)\n",
        "// General form (base):\n",
        "//   Um = Sum_j[mu_j(t,rho_vac[SCm])/r_j*(1-exp(-gamma*t)*cos(pi*t_n))*phi^j]\n",
        "//       * P_SCm * E_react * (1+1e13*f_Heaviside) * (1+f_quasi) * exp(-[SSq])\n",
        "//       * TYPE_SPECIFIC_FACTOR\n",
        "// ======================================================\n",
        "std::vector<BuoyancyEquation> buildSectionF() {\n",
        "    return {\n",
    ]
    for entry in _F_DATA:
        id_, name, expr, ctx, note, num, units, system = entry
        lines.append(
            "        { "
            + str(id_) + ", "
            + '"' + _esc(name)   + '", '
            + '"' + _esc(expr)   + '", '
            + '"' + _esc(ctx)    + '", '
            + '"' + _esc(note)   + '", '
            + repr(float(num))   + ", "
            + '"' + _esc(units)  + '", '
            + '"' + _esc(system) + '", "F" },\n'
        )
    lines.append("    };\n}\n")
    return "".join(lines)
