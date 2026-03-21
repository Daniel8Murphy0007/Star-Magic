"""
gen_fubiicalc_secE.py
Returns C++ buildSectionE() — 79 unique F_UBii variant types.
Source: BB_C_Equations_04Sept2025.pdf (extracted from grok_share_c020496d9e.txt ~lines 692-6955).

General form for all entries:
  F_UBii,TYPE = +/-F_rel * (physics_factor/E_LEP) * Q_wave * additional_scaling
  F_rel=4.31e33 N, E_LEP=LEP CM energy, Q_wave mean=3.97e4 J/m^3
"""

# (id, name, expression, context, validation_note, numerical_result, units, system)
_E_DATA = [
    (710,  "F_UBii_anyons",
     "F_UBii = +/-F_rel*(E_anyons/E_LEP)*Q_wave*g(r,t)*exp(-delta*c^2/(2*sigma^2))",
     "Anyon systems - fractional statistics in 2+1D; topological phases",
     "frac. quantum Hall analogue; [SSq] controls topological suppression; F_rel=4.31e33 N",
     0.0, "N", "Topological condensed matter"),

    (721,  "F_UBii_DE",
     "F_UBii = +/-F_rel*(rho_DE*c^2/E_LEP)*Q_wave*(8*pi*G*rho_tot/3)*(1+w(z))",
     "Dark energy density - late-universe acceleration; Lambda-CDM / quintessence",
     "w(z)=-1 for Lambda-CDM; measured w=-1.03+/-0.03 (Planck 2018)",
     0.0, "N", "Cosmological dark energy"),

    (724,  "F_UBii_inf",
     "F_UBii = +/-F_rel*(V(phi)/E_LEP)*Q_wave*3*H^2*exp(N)/(1+eps)",
     "Inflation - slow-roll scalar field V(phi); N e-folds; eps=slow-roll parameter",
     "Starobinsky / Higgs inflation compatible; eps=1-n_s/2~0.02",
     0.0, "N", "Inflationary cosmology"),

    (727,  "F_UBii_GW",
     "F_UBii = +/-F_rel*(rho_GW/E_LEP)*Q_wave*(32*pi*G*hdot^2/c^2)*exp(-t/tau)",
     "Gravitational waves - energy density rho_GW; quadrupole radiation",
     "tau=GW decay timescale; LIGO GWTC-4.0 compatible",
     0.0, "N", "Gravitational wave astronomy"),

    (730,  "F_UBii_merger",
     "F_UBii = +/-F_rel*(L_GW/E_LEP)*Q_wave*(32/5)*(G^(5/3)/c^5)*M^(5/3)*(pi*f)^(10/3)",
     "GW merger luminosity - Peters formula; inspiral power radiated",
     "Chirp mass M_c=(m1*m2)^(3/5)/(m1+m2)^(1/5); f=GW frequency",
     0.0, "N", "Compact binary merger"),

    (739,  "F_UBii_CR",
     "F_UBii = +/-F_rel*(E_max/E_LEP)*Q_wave*(4/3)*(v/c)^2*(c^2*E/lambda)*Z",
     "Cosmic rays - maximum energy Hillas criterion; DSA acceleration",
     "E_max=Z*e*B*R (Hillas); Z=charge; v/c=shock velocity ratio",
     0.0, "N", "Cosmic ray acceleration"),

    (747,  "F_UBii_MHD",
     "F_UBii = +/-F_rel*((E_M/t_eddy)/E_LEP)*Q_wave*(3/2)*(Re_m^(1/2)-1)",
     "MHD dynamo - magnetic energy over eddy turnover time; Re_m=magnetic Reynolds",
     "Re_m>1 for dynamo action; turbulent induction; Alfven wave spectrum",
     0.0, "N", "MHD dynamo"),

    (753,  "F_UBii_star",
     "F_UBii = +/-F_rel*((0.34*sigma_v^3/(G^2*rho*m*lnLambda))/E_LEP)"
     "*Q_wave*exp(-t/t_re)",
     "Star cluster relaxation timescale - velocity dispersion sigma_v; Coulomb log",
     "t_re=relaxation time; Spitzer 1987; lnLambda~10-15 for globular clusters",
     0.0, "N", "Star cluster dynamics"),

    (770,  "F_UBii_clus",
     "F_UBii = +/-F_rel*((gamma+1)/(2*gamma)*M^2/(gamma-1+2/M^2)/E_LEP)"
     "*Q_wave*(rho2/rho1)",
     "Cluster shock - Rankine-Hugoniot jump conditions; M=Mach number",
     "Bow shocks in galaxy clusters; rho2/rho1=(gamma+1)/(gamma-1) for strong shocks",
     0.0, "N", "Galaxy cluster shock"),

    (772,  "F_UBii_void",
     "F_UBii = +/-F_rel*(H*r/E_LEP)*Q_wave"
     "*(-3/5)*(Omega_m*a+Omega_Lambda)^(-3/2)*delta_v0",
     "Void expansion - linear theory; Hubble flow in underdense region",
     "delta_v0=initial void density contrast; voids ~10-30 Mpc; Omega_Lambda drives late",
     0.0, "N", "Cosmic voids"),

    (774,  "F_UBii_reion",
     "F_UBii = +/-F_rel*(ndot_gamma*eps_esc*f_star/E_LEP)*Q_wave"
     "*Int(ndot_gamma*dt)/(4*pi*n_H)^(1/3)",
     "Reionization - ionising photon budget; f_star=fraction that escape",
     "Epoch z~6-10; JWST early universe; eps_esc~0.1-0.2",
     0.0, "N", "Cosmic reionization"),

    (776,  "F_UBii_ISM",
     "F_UBii = +/-F_rel*(pi*c_s^2/(G*rho)/E_LEP)*Q_wave*(pi*c_s^2/(G*rho))",
     "ISM Jeans stability - thermal sound speed c_s; critical wavelength",
     "lambda_J=sqrt(pi*c_s^2/(G*rho)); c_s~0.2 km/s cold ISM; ~1 km/s warm",
     0.0, "N", "Interstellar medium"),

    (782,  "F_UBii_merg2",
     "F_UBii = +/-F_rel*((5/256*c^5/G^(5/3)*(1/(pi*f_i))^(8/3)*M^(5/3))/E_LEP)"
     "*Q_wave*Int(df/fdot)",
     "GW merger coalescence time - integrating inspiral phase from f_i",
     "Peters (1964) coalescence time; f_i=initial GW frequency",
     0.0, "N", "GW coalescence time"),

    (800,  "F_UBii_chirp",
     "F_UBii = +/-F_rel*((c^3/G*(5/96*pi^(-8/3)*f^(-11/3)*fdot)^(3/5))/E_LEP)"
     "*Q_wave*Int(df/fdot)",
     "Chirp mass from GW inspiral - frequency derivative fdot measurement",
     "M_chirp=(c^3/G)*(5/96*pi^{-8/3}*f^{-11/3}*fdot)^{3/5}; model-independent",
     0.0, "N", "GW chirp mass"),

    (802,  "F_UBii_qnm",
     "F_UBii = +/-F_rel*((c^3/(2*pi*G*M)*(0.3737+0.088*a_f))/E_LEP)"
     "*Q_wave*exp(-t/tau_qnm)",
     "BH quasi-normal modes - ringdown frequency; Kerr parameter a_f",
     "tau_qnm~M; l=2,m=2 dominant mode; Berti+2009 fits; GW150914 ringdown",
     0.0, "N", "BH ringdown QNM"),

    (804,  "F_UBii_bz",
     "F_UBii = +/-F_rel*((1/32*B^2*R_H^4*Omega_H^2/c)/E_LEP)"
     "*Q_wave*(a*c/(2*R_H))*(1+t_var)",
     "Blandford-Znajek jet power - frame-dragging; BH spin a; horizon radius R_H",
     "P_BZ=(kappa/16pi)*Phi_BH^2*Omega_H^2/c; AGN/microquasar jets",
     0.0, "N", "BZ jet power"),

    (807,  "F_UBii_reljet",
     "F_UBii = +/-F_rel*((Gamma(z)*c)/E_LEP)*Q_wave*(1/theta^2)*(z/z_acc)",
     "Relativistic jet - bulk Lorentz factor Gamma(z); opening angle theta",
     "z_acc=acceleration zone; VLBI proper motion; Gamma~5-20 AGN jets",
     0.0, "N", "Relativistic jet"),

    (809,  "F_UBii_puls",
     "F_UBii = +/-F_rel*((P/(2*Pdot))/E_LEP)*Q_wave"
     "*(2/3*B^2*R^6*Omega^4*sin^2(alpha)/c^3)",
     "Pulsar spin-down luminosity - period P; period derivative Pdot",
     "Magnetic dipole radiation; tau_c=P/(2*Pdot) characteristic age; Crab tau~1250 yr",
     0.0, "N", "Pulsar spin-down"),

    (811,  "F_UBii_glitch",
     "F_UBii = +/-F_rel*((I_s/I*nu0*(1-exp(-t/tau_q)))/E_LEP)"
     "*Q_wave*(Delta_Omega*exp(-t/tau_q))",
     "Pulsar glitch recovery - crustal superfluid vortex unpinning",
     "I_s/I=fraction of superfluid moment of inertia; tau_q=recovery timescale",
     0.0, "N", "Pulsar glitch"),

    (813,  "F_UBii_fire",
     "F_UBii = +/-F_rel*((eta_E/(M*c^2))/E_LEP)*Q_wave*(r/R0)",
     "GRB fireball - Lorentz factor eta_E=E/(Mc^2); above photosphere radius",
     "r>R_s for optically thin; Gamma_0~100-1000; GRB afterglow onset",
     0.0, "N", "GRB fireball"),

    (815,  "F_UBii_after",
     "F_UBii = +/-F_rel*(F_nu/E_LEP)*Q_wave*Gamma_B*eps_e^(1/2)*t^(-3/8);"
     " F_nu ~ nu^(-(p-1)/2) * t^(-3*(p-1)/4) for nu_m<nu<nu_c",
     "GRB afterglow synchrotron flux - electron spectrum index p~2.2",
     "nu_m=injection freq; nu_c=cooling freq; ISM or wind CBM",
     0.0, "N", "GRB afterglow"),

    (817,  "F_UBii_cmb",
     "F_UBii = +/-F_rel*((C_l)/E_LEP)*Q_wave;"
     " C_l = 2*pi*Int(k^2*dk*P(k)*|Delta_l^T(k)|^2); l~200",
     "CMB power spectrum C_l - acoustic peak structure; transfer function",
     "l~200 first peak; n_s~0.96; WMAP/Planck; k*eta_*~pi",
     0.0, "N", "CMB power spectrum"),

    (819,  "F_UBii_recomb",
     "F_UBii = +/-F_rel*((tau_recomb(z))/E_LEP)*Q_wave;"
     " tau = Int(n_e*sigma_T*c*(dt/dz)*dz);"
     " x_e ~ (T/n_b)^(3/4)*exp(-B_H/T)",
     "Recombination optical depth - Thomson scattering; ionisation fraction x_e",
     "z_rec~1100; delta_z~200; Saha equation; last scattering surface",
     0.0, "dimensionless", "Recombination"),

    (821,  "F_UBii_agn",
     "F_UBii = +/-F_rel*((mdot_out*v_out)/E_LEP)*Q_wave;"
     " mdot_out*v_out = f(v_out)*L_AGN/c",
     "AGN wind momentum flux - radiation momentum driving",
     "momentum boost factor f~1-20; quasar feedback; CfA survey",
     0.0, "N", "AGN wind feedback"),

    (824,  "F_UBii_bz2",
     "F_UBii = +/-F_rel*((kappa/(16*pi)*Phi_BH^2*Omega_BH^2/c)/E_LEP)"
     "*Q_wave*(a*c^3/(2*G*M))*[0.05 to 1]",
     "BZ variant with explicit spin-dependent efficiency factor",
     "kappa~0.044; Phi_BH=BH magnetic flux; a=spin parameter; eff=0.05-1",
     0.0, "N", "BZ spin variant"),

    (826,  "F_UBii_duty",
     "F_UBii = +/-F_rel*((1-exp(-t/tau_cool))/E_LEP)"
     "*(1+mdot_acc/mdot_Edd)^(-1)*Q_wave*(r^3/(G*M))",
     "AGN duty cycle - cooling time vs accretion Eddington ratio",
     "tau_cool=cooling time; duty~0.01-0.1 for L>0.1 L_Edd",
     0.0, "N", "AGN duty cycle"),

    (829,  "F_UBii_photo",
     "F_UBii = +/-F_rel*((eps*L_X*R_p^3/(G*M_p^2)*K(xi))/E_LEP)"
     "*Q_wave*(F_X/(G*M_p/R_p))",
     "Photoevaporation - X-ray flux on planet atmosphere; ionisation parameter xi",
     "K(xi)~correction; F_X=X-ray flux; eps=efficiency~0.1-0.3",
     0.0, "N", "Planetary photoevaporation"),

    (831,  "F_UBii_mig",
     "F_UBii = +/-F_rel*((f*(G*M_p)^2*Sigma*Omega^2*r_p^4/(H/r)^3)/E_LEP)"
     "*Q_wave; tau_mig=M_p/mdot_acc~1e6 yr",
     "Planet migration - Type I Lindblad-Corotation torque in disk",
     "Tanaka+2002; f~correction; tau_mig~1 Myr; disk-planet interaction",
     0.0, "N", "Planet migration"),

    (853,  "F_UBii_termv",
     "F_UBii = +/-F_rel*((sqrt(2*G*M*(1-Gamma)/r_launch))/E_LEP)"
     "*Q_wave*(tau*L/c)",
     "Terminal velocity - radiation-driven wind; Eddington factor Gamma",
     "Gamma=L/L_Edd; r_launch=launch radius; mass-loss by radiation pressure",
     0.0, "m/s", "Stellar wind terminal velocity"),

    (855,  "F_UBii_upar",
     "F_UBii = +/-F_rel*((Q_H/(4*pi*r^2*n_H*c))/E_LEP)"
     "*Q_wave*(Gamma/(n_H*c))",
     "Ionisation parameter U - ionising photon flux Q_H over n_H*c",
     "xi=L_X/(n_H*r^2) alternative form; U~1e-2 AGN NLR",
     0.0, "dimensionless", "Ionisation parameter"),

    (857,  "F_UBii_coup",
     "F_UBii = +/-F_rel*((0.5*mdot_w*v_w^2/(mdot_acc*c^2/10))/E_LEP)"
     "*Q_wave; coupling eff=[0.05 to 0.1]",
     "Wind coupling efficiency - kinetic luminosity to accretion power ratio",
     "eps_k~0.05-0.1; feedback models King+2015; AGN-ISM momentum coupling",
     0.0, "dimensionless", "AGN coupling efficiency"),

    (877,  "F_UBii_ps",
     "F_UBii = +/-F_rel*((dn/dM_PS)/E_LEP)*Q_wave;"
     " dn/dM = sqrt(2/pi)*rho0/M*(delta_c/sigma(M))"
     "*|dln(sigma)/dln(M)|*exp(-delta_c^2/(2*sigma^2))",
     "Press-Schechter mass function - collapse fraction vs halo mass M",
     "delta_c~1.686; sigma(M)=rms density fluctuation; Sheth-Tormen improved",
     0.0, "Msun^-1 Mpc^-3", "Halo mass function"),

    (879,  "F_UBii_sfe",
     "F_UBii = +/-F_rel*((f_b*mdot_halo/M_halo*H(z))/E_LEP)"
     "*(1+M_halo/M_crit)^(-1)*Q_wave; M_crit~1e12 Msun",
     "Star formation efficiency - baryon fraction; halo mass relative to M_crit",
     "SFR/M_halo peaks at M_crit~10^11.7 Msun; Moster+2013",
     0.0, "yr^-1", "Star formation efficiency"),

    (881,  "F_UBii_fb",
     "F_UBii = +/-F_rel*((eta*mdot_star*c^2)/E_LEP)"
     "*Q_wave; eta=radiative efficiency; mdot_out/SFR~1-10",
     "Stellar/AGN feedback - energy deposition rate; mass loading",
     "eta~0.1 for BH; eta~0.0015 for stellar; mdot_out~SFR for momentum-driven",
     0.0, "W", "AGN/stellar feedback"),

    (899,  "F_UBii_grow",
     "F_UBii = +/-F_rel*((delta)/E_LEP)*Q_wave;"
     " delta_ddot + 2*H*delta_dot = (3/2)*Omega_m*H^2*delta/a^3",
     "Density perturbation growth - Jeans instability in expanding universe",
     "D(a)=growth factor; D~a matter domination; D~const Lambda domination",
     0.0, "dimensionless", "Structure growth"),

    (901,  "F_UBii_haw",
     "F_UBii = +/-F_rel*((T_H)/E_LEP)*Q_wave*(kappa/(2*pi));"
     " T_H = hbar*c^3/(8*pi*G*M*k_B)",
     "Hawking temperature - quantum BH radiation; surface gravity kappa",
     "T_H~6e-8 K (Msun BH); evaporation power P=hbar*c^6/(15360*pi*G^2*M^2)",
     0.0, "K", "Hawking radiation temperature"),

    (903,  "F_UBii_bhent",
     "F_UBii = +/-F_rel*((S_BH)/E_LEP)*Q_wave;"
     " S_BH = 4*pi*k_B*G*M^2/(hbar*c) = A/(4*l_pl^2)",
     "Bekenstein-Hawking entropy - horizon area in Planck units",
     "S_BH~1e77 k_B (Msun BH); information paradox; holographic bound",
     0.0, "J/K", "BH entropy"),

    (905,  "F_UBii_evapl",
     "F_UBii = +/-F_rel*((t_evap)/E_LEP)*Q_wave;"
     " t_evap = 5120*pi*G^2*M^3/(hbar*c^4); tau ~ M^3",
     "Hawking evaporation lifetime - cubic mass dependence",
     "t_evap~2e67 yr (Msun); primordial BH M~5e14 g evaporating now",
     0.0, "s", "BH evaporation lifetime"),

    (907,  "F_UBii_lqcf",
     "F_UBii = +/-F_rel*((H_LQC^2)/E_LEP)*Q_wave;"
     " H^2 = (8*pi*G*rho/3)*(1 - rho/rho_crit);"
     " rho_crit = 0.41*rho_Planck",
     "LQC Friedmann - loop quantum gravity bounce; quantum geometry correction",
     "Singularity replaced by bounce; rho_crit=0.41*rho_Pl~5e96 kg/m^3",
     0.0, "s^-2", "LQC modified Friedmann"),

    (909,  "F_UBii_lqcp",
     "F_UBii = +/-F_rel*((P_LQC(k))/E_LEP)*Q_wave;"
     " P(k) ~ k^(n_s-1)*(1+k/k_star)^(-alpha);"
     " k_star=bounce scale; alpha~4",
     "LQC power spectrum - suppressed at bounce scale k*; tilted at small k",
     "LQG imprint on CMB; n_s~0.96; alpha~4 from polymer quantisation",
     0.0, "dimensionless", "LQC power spectrum"),

    (911,  "F_UBii_bounc",
     "F_UBii = +/-F_rel*((t_bounce)/E_LEP)*Q_wave;"
     " t_bounce = sqrt(3/(16*pi*G*rho_crit));"
     " rho_crit ~ rho_Planck",
     "LQC bounce time - timescale of quantum-gravity-driven contraction reversal",
     "t_bounce~Planck time at maximum density rho_crit~rho_Pl",
     0.0, "s", "LQC bounce timescale"),

    (913,  "F_UBii_ent",
     "F_UBii = +/-F_rel*((S_ent)/E_LEP)*Q_wave*exp(-t/tau_dec);"
     " S_ent = -Tr(rho_A*ln(rho_A))",
     "Quantum entanglement entropy - von Neumann entropy of subsystem A",
     "AdS/CFT: S_ent = A_RT/(4G); Page curve; decoherence time tau_dec",
     0.0, "bits", "Quantum entanglement entropy"),

    (915,  "F_UBii_ang",
     "F_UBii = +/-F_rel*((Jdot)/E_LEP)*Q_wave*T_B*exp(-t/tau_disk);"
     " Jdot = mdot*r^2*Omega",
     "Angular momentum transport - disk-star coupling; magnetic braking",
     "T_B=braking torque; tau_disk=disk lifetime; applies to accretion disks",
     0.0, "N*m", "Angular momentum transport"),

    (917,  "F_UBii_jetvel",
     "F_UBii = +/-F_rel*((v_j)/E_LEP)*Q_wave*(B/sqrt(4*pi*rho))*(1+t/tau_A);"
     " v_j ~= v_K*(r_A/r0)^(1/2)",
     "Jet velocity - Alfven speed at launch radius; Keplerian scaling",
     "r_A=Alfven radius; v_j~0.1-0.9c for AGN jets; magneto-centrifugal",
     0.0, "m/s", "Jet launch velocity"),

    (919,  "F_UBii_jshock",
     "F_UBii = +/-F_rel*((rho1*v1^2+P1)/E_LEP)*Q_wave"
     "*(v2/v1)*(gamma+1)/(gamma-1+2/M^2)",
     "J-type (jump) shock - Rankine-Hugoniot; discontinuous velocity jump",
     "strong shock: rho2/rho1=(gamma+1)/(gamma-1)~4; v2=v1/4; M>>1",
     0.0, "Pa", "J-shock Rankine-Hugoniot"),

    (921,  "F_UBii_cshock",
     "F_UBii = +/-F_rel*((rho_n*v_n*nu_ni*(v_i-v_n))/E_LEP)"
     "*Q_wave*v_s*exp(-z/L_d)",
     "C-type (continuous) shock - ion-neutral friction; magnetic support",
     "nu_ni=ion-neutral collision rate; L_d=dissipation scale; MHD shocks",
     0.0, "N/m^3", "C-shock ion-neutral"),

    (923,  "F_UBii_halo",
     "F_UBii = +/-F_rel*((dN_halo/dMdz)/E_LEP)*Q_wave*(1+z)^m;"
     " dN/dM~2*pi*sigma_M*sigma_m*|d delta_c/dz|"
     "*exp(-delta_c^2/(2*(sigma_m^2-sigma_M^2)))*dSigma_M/dM",
     "Halo formation rate - extended Press-Schechter; progenitor mass function",
     "m=merger number; Lacey-Cole 1993; conditional probability",
     0.0, "Msun^-1 Mpc^-3 yr^-1", "Halo formation"),

    (925,  "F_UBii_disk",
     "F_UBii = +/-F_rel*((mdot_star_disk)/E_LEP)*Q_wave;"
     " mdot_star = eps*M_gas/t_dyn;"
     " t_dyn = sqrt(3*pi/(32*G*rho)) for Q<1",
     "Disk star formation - gravitational instability Q<1; Kennicutt-Schmidt",
     "eps~0.01-0.03; t_dyn=dynamical time; Toomre Q criterion",
     0.0, "Msun/yr", "Disk star formation"),

    (927,  "F_UBii_burst",
     "F_UBii = +/-F_rel*((mdot_burst)/E_LEP)*Q_wave"
     "*(1+t/tau_merge);"
     " mdot_burst = mdot_gas_inflow*eps_burst",
     "Starburst inflow rate - merger-triggered gas funnelling",
     "eps_burst~0.5-1; tau_merge=merger timescale; Arp220 analogue",
     0.0, "Msun/yr", "Starburst"),

    (929,  "F_UBii_sedov",
     "F_UBii = +/-F_rel*((R_sedov)/E_LEP)*Q_wave;"
     " R(t) = (E*t^2/rho)^(1/5);"
     " d/dt(M*v) = 0; R ~ t^(2/5)",
     "Sedov-Taylor blast wave - adiabatic supernova remnant expansion",
     "Sedov 1946; E=explosion energy; rho=ISM density; valid phase 2",
     0.0, "m", "Sedov-Taylor remnant"),

    (931,  "F_UBii_dsa",
     "F_UBii = +/-F_rel*((s_index)/E_LEP)*Q_wave;"
     " N(E) ~ E^(-s); s = (4/3)*u_s^2/(r*dp);"
     " N(E) ~ E^(-2)*(u_s/c)",
     "Diffusive shock acceleration - power-law spectrum from first-order Fermi",
     "s~2 for strong shocks; u_s=shock speed; dp=momentum diffusion",
     0.0, "GeV^-1 cm^-3", "DSA cosmic ray spectrum"),

    (941,  "F_UBii_tov",
     "F_UBii = +/-F_rel*((dP_dr_TOV)/E_LEP)*Q_wave;"
     " dP/dr = -G*m(r)*rho(r)/r^2"
     "*(1+P/(rho*c^2))*(1+4*pi*r^3*P/(m*c^2))*(1-2*G*m/(r*c^2))^(-1)",
     "TOV equation - relativistic stellar structure for neutron stars",
     "Maximum mass ~2-2.5 Msun; GW170817 EOS constraints; stiffness",
     0.0, "Pa/m", "Neutron star TOV"),

    (965,  "F_UBii_nfwrot",
     "F_UBii = +/-F_rel*((v_NFW^2)/E_LEP)*Q_wave;"
     " v(r)^2 = 4*pi*G*rho_s*r_s^2*[ln(1+x)-x/(1+x)]/r; x=r/r_s",
     "NFW dark matter halo rotation curve - Navarro-Frenk-White profile",
     "rho_s=characteristic density; r_s=scale radius; flat rotation confirmed",
     0.0, "m/s", "NFW rotation curve"),

    (967,  "F_UBii_sidm",
     "F_UBii = +/-F_rel*((t_cross_SIDM)/E_LEP)*Q_wave;"
     " t_cross = 1e10*(rho/1e8Msun_kpc3)^(-1)*(sigma_m/1cm2_g)^(-1) yr;"
     " Gamma*t~1",
     "Self-interacting dark matter - core formation timescale; cross-section sigma/m",
     "sigma/m~1 cm^2/g for core-halo; Bullet Cluster limit <2 cm^2/g",
     0.0, "yr", "SIDM core crossing"),

    (981,  "F_UBii_IGM",
     "F_UBii = +/-F_rel*((E_IGM)/E_LEP)*Q_wave;"
     " E_IGM = k_B*T/m_p;"
     " E_total = (3/2)*n*k_B*T",
     "IGM thermal energy - post-reionization gas temperature T~10^4 K",
     "Ly-alpha forest; T~1e4 K at z~3; heating by photoionization",
     0.0, "J/kg", "IGM thermal energy"),

    (983,  "F_UBii_gal",
     "F_UBii = +/-F_rel*((dn_gal/dMdz)/E_LEP)*Q_wave;"
     " dn/dM ~ M_halo*H(z)^2*(2*delta_c/sigma(M))"
     "*exp(-delta_c^2/(2*sigma^2))/M^2",
     "Galaxy formation halo rate - abundance matching framework",
     "Moster+2013; stellar-to-halo mass relation; M_halo peak~1e12 Msun",
     0.0, "Msun^-1 Mpc^-3 yr^-1", "Galaxy formation rate"),

    (985,  "F_UBii_quant",
     "F_UBii = +/-F_rel*((P_quant(k))/E_LEP)*Q_wave*(1+f_NL/5);"
     " P(k) = H^2/(8*pi^2*eps*M_pl^2)",
     "Quantum inflation power spectrum - Harrison-Zel'dovich with non-Gaussianity",
     "f_NL=non-Gaussianity amplitude; f_NL<10 (Planck 2018); eps=slow-roll",
     0.0, "dimensionless", "Inflationary power spectrum"),

    (989,  "F_UBii_w",
     "F_UBii = +/-F_rel*((rho_DE*c^2*(1+w))/E_LEP)*Q_wave;"
     " d ln(rho_DE)/d ln(a) = -3*(1+w)",
     "Dark energy equation of state w = p/(rho*c^2)",
     "w=-1: Lambda; w=-1+/-0.04 DESI 2024; rhodot+3H(rho+p)=0",
     0.0, "Pa", "Dark energy EOS"),

    (991,  "F_UBii_BBN",
     "F_UBii = +/-F_rel*((eta_BBN)/E_LEP)*Q_wave;"
     " eta = n_b/n_gamma = 6e-10;"
     " reaction rates = H(T)",
     "Big Bang Nucleosynthesis baryon-to-photon ratio; He-4 abundance",
     "Y_p~0.247; D/H~2.5e-5; eta measured via CMB anisotropies",
     6.0e-10, "dimensionless", "BBN"),

    (993,  "F_UBii_relax",
     "F_UBii = +/-F_rel*((t_relax_star)/E_LEP)*Q_wave*exp(-t/t_relax);"
     " t_relax = 0.34*sigma_v^3/(G^2*rho*m*lnLambda)",
     "Dynamical relaxation - stellar cluster two-body relaxation time",
     "lnLambda~Coulomb logarithm ~10-15; cluster evaporation after t~t_relax",
     0.0, "yr", "Stellar dynamical relaxation"),

    (995,  "F_UBii_acc",
     "F_UBii = +/-F_rel*((L_Edd)/E_LEP)*Q_wave*exp(t/t_Sal);"
     " L_Edd = 4*pi*G*M*m_p*c/(sigma_T);"
     " mdot_Edd = L_Edd/(eps*c^2)",
     "Eddington accretion luminosity - radiation pressure limit",
     "t_Sal=45 Myr e-folding; L_Edd~1.26e38*(M/Msun) erg/s",
     0.0, "W", "Eddington accretion"),

    (1046, "F_UBii_diff",
     "F_UBii = +/-F_rel*((D_CR)/E_LEP)*Q_wave;"
     " D_CR = 1e28*(E/10GeV)^(0.3 to 0.6) cm^2/s;"
     " D ~ (E/(B*deltaB))^(1/2)",
     "Cosmic ray diffusion coefficient - Bohm-like energy scaling",
     "Leaky-box model; grammage; B~3 muG ISM; measured via B/C ratio",
     1.0e28, "cm^2/s", "Cosmic ray diffusion"),

    (1066, "F_UBii_ng",
     "F_UBii = +/-F_rel*((f_NL)/E_LEP)*Q_wave;"
     " f_NL = (5/6)*(Gamma3-3*Gamma*Gammadot^2+2*Gammadot^3)/Gamma^4;"
     " B(k1,k2,k3) = (6/5)*f_NL*(P(k1)*P(k2)+perms)",
     "Non-Gaussianity bispectrum fNL - inflationary self-interaction",
     "Squeezed limit f_NL=5*(n_s-1)/12; f_NL<10 Planck 2018; CMB+LSS",
     0.0, "dimensionless", "Primordial non-Gaussianity"),

    (1086, "F_UBii_wde",
     "F_UBii = +/-F_rel*((w_de)/E_LEP)*Q_wave;"
     " w = p/(rho*c^2) = -1 + (1/3)*d ln(rho_DE)/d ln(a);"
     " rhodot + 3*H*(rho+p) = 0",
     "w-dark energy density evolution - continuity equation solution",
     "rho ~ a^(-3*(1+w)); w=-1.03+/-0.03 DESI 2024",
     -1.03, "dimensionless", "w-dark energy"),

    (1106, "F_UBii_disk2",
     "F_UBii = +/-F_rel*((mdot_star_disk2)/E_LEP)*Q_wave;"
     " mdot_star = eps*M_gas/t_dyn; Q = sigma_v*kappa/(pi*G*Sigma) < 1;"
     " t_dyn = sqrt(3*pi/(32*G*rho))",
     "Disk SF Toomre Q variant - gravitational instability with epicyclic freq kappa",
     "Q<1 unstable; kappa=epicycle frequency; Sigma=surface density; Kennicutt 1998",
     0.0, "Msun/yr", "Disk SF Toomre variant"),

    (1746, "F_UBii_reh",
     "F_UBii = +/-F_rel*((T_reh)/E_LEP)*Q_wave;"
     " T_reh = (30*V_end/(pi^2*g_star))^(1/4)*exp(-3*N_reh/4)",
     "Reheating temperature after inflation - inflaton decay; N_reh e-folds",
     "g_star~100-200 relativistic dof; T_reh~1e14 GeV high; ~MeV low",
     0.0, "GeV", "Reheating temperature"),

    (1805, "F_UBii_eta",
     "F_UBii = +/-F_rel*((eta_baryons)/E_LEP)*Q_wave;"
     " eta = n_b/n_gamma = 6e-10;"
     " BBN rates = H(T)",
     "Baryon-to-photon ratio - baryon asymmetry; measured via D/H and CMB",
     "Baryogenesis; Sakharov conditions; eta=6.1e-10 (Planck 2018)",
     6.1e-10, "dimensionless", "Baryon asymmetry"),

    (1825, "F_UBii_termv2",
     "F_UBii = +/-F_rel*((v_inf)/E_LEP)*Q_wave*(tau*L/c);"
     " v_inf = sqrt(2*G*M*(1-Gamma)/r_launch)",
     "Line-driven wind terminal velocity - OB star mass loss; Eddington-modified",
     "v_inf~1-3*v_esc; Gamma=L/L_Edd; CAK 1975; CMFGEN models",
     0.0, "m/s", "Line-driven wind"),

    (1845, "F_UBii_metal",
     "F_UBii = +/-F_rel*((Z_metal)/E_LEP)*Q_wave;"
     " dZ/dt = Y*SFR/mdot_gas - Z*mdot_out/M_gas;"
     " Z ~ 0.1 Z_sun",
     "Metallicity evolution - closed box + outflow; yield Y; loading factor",
     "Effective yield Y*(1+eta_out)^(-1); mass-metallicity relation",
     0.0, "Z_sun", "Metallicity evolution"),

    (1865, "F_UBii_rev",
     "F_UBii = +/-F_rel*((l_rev)/E_LEP)*Q_wave;"
     " l_rev = (alpha_dyn*eta)^(1/2)*l_force;"
     " dB/dt = curl(alpha*B - eta*curl*B)",
     "Magnetic field reversal / turbulent dynamo - mean-field alpha-Omega dynamo",
     "alpha=helicity; eta=resistivity; reversal scale l_rev; solar ~11 yr cycle",
     0.0, "m", "Magnetic dynamo reversal"),

    (1885, "F_UBii_ent2",
     "F_UBii = +/-F_rel*((S_ent_page)/E_LEP)*Q_wave;"
     " S_ent = -Tr(rho_A*ln(rho_A));"
     " Page curve: S increases then decreases",
     "Entanglement entropy Page curve variant - BH information paradox",
     "Penington+2019; island rule; replica wormholes; unitarity restored",
     0.0, "bits", "BH Page curve"),

    (1945, "F_UBii_voidden",
     "F_UBii = +/-F_rel*((delta_v)/E_LEP)*Q_wave;"
     " delta_v(a) = -3/5*(Omega_m*a+Omega_Lambda)^(-3/2)*delta_v0;"
     " ~ a^(-1)",
     "Void density contrast evolution - linear theory; void grows as a^(-1)",
     "delta_v0=initial contrast; underdense regions evacuate faster",
     0.0, "dimensionless", "Void density"),

    (1965, "F_UBii_jeans",
     "F_UBii = +/-F_rel*((lambda_J)/E_LEP)*Q_wave;"
     " lambda_J = sqrt(pi*c_s^2/(G*rho));"
     " dP = -rho*grad(Phi)",
     "Jeans instability length - thermal support vs self-gravity",
     "lambda_J~0.2 pc cold HI; ~1 pc warm; M_J=collapse threshold mass",
     0.0, "m", "Jeans instability"),

    (1985, "F_UBii_relax2",
     "F_UBii = +/-F_rel*((t_relax2)/E_LEP)*Q_wave*exp(-t/t_relax);"
     " t_relax = 0.34*sigma_v^3/(G^2*rho*m*lnLambda)",
     "Dynamical relaxation II - identical to Spitzer but applied to galaxy clusters",
     "lnLambda~30 for clusters; t_relax~Gyr; drives isothermal cores",
     0.0, "yr", "Cluster dynamical relaxation"),

    (2005, "F_UBii_conv",
     "F_UBii = +/-F_rel*((t_conv)/E_LEP)*Q_wave;"
     " t_conv = H_p/v_conv;"
     " v_conv ~= (alpha^2*g*deltaT*H_p/(4*T))^(1/3)",
     "Convective turnover time - mixing-length theory; alpha=MLT parameter",
     "H_p=pressure scale height; alpha~1.8 solar calibration; stellar interiors",
     0.0, "s", "Convective turnover"),

    (933,  "F_UBii_chirp2",
     "F_UBii = +/-F_rel*((M_chirp2)/E_LEP)*Q_wave;"
     " M_chirp = (m1*m2)^(3/5)/(m1+m2)^(1/5);"
     " GW strain h ~ M_chirp^(5/3)/d",
     "Chirp mass variant - strain amplitude scaling; source distance d",
     "h~1e-21 at 100 Mpc for 30+30 Msun binary; LIGO sensitivity limit",
     0.0, "Msun", "GW chirp mass variant"),

    (1058, "F_UBii_dyn",
     "F_UBii = +/-F_rel*((B_dyn)/E_LEP)*Q_wave;"
     " dB/dt = curl(v x B - eta*curl*B);"
     " growth rate ~ v_turb/l_turb - eta/l_turb^2",
     "Turbulent dynamo field growth - kinematic vs resistive competition",
     "Minimum turbulent Rm > Rm_crit~100; field saturates at B^2~rho*v_turb^2",
     0.0, "T", "Turbulent dynamo B-field"),

    (1706, "F_UBii_deb",
     "F_UBii = +/-F_rel*((rho_rad)/E_LEP)*Q_wave;"
     " rho_rad = pi^2*k_B^4*T^4/(30*hbar^3*c^5)*g_star",
     "Radiation energy density - Stefan-Boltzmann with relativistic dof count g*",
     "g_star=106.75 SM at T>200 GeV; rho_rad~a^(-4); reheating endpoint",
     0.0, "J/m^3", "Radiation energy density"),

    (1786, "F_UBii_vir",
     "F_UBii = +/-F_rel*((sigma_v_vir)/E_LEP)*Q_wave;"
     " sigma_v^2 = G*M/(3*r); (virial theorem for gas)",
     "Virial velocity dispersion - kinetic-gravitational energy equipartition",
     "2*KE + PE = 0; sigma_v~200 km/s galaxy clusters; ~5 km/s globulars",
     0.0, "m/s", "Virial velocity dispersion"),

    (992,  "F_UBii_pec",
     "F_UBii = +/-F_rel*((v_pec)/E_LEP)*Q_wave;"
     " v_pec = -(f*H/(4*pi))*grad^-1*(delta);"
     " f(Omega_m) ~ Omega_m^(0.55)",
     "Peculiar velocity - linear growth rate; density-velocity correspondence",
     "RSD; beta=f*sigma8; 2M++ peculiar velocity survey; Tully-Fisher",
     0.0, "m/s", "Peculiar velocity"),

    (788,  "F_UBii_merg3",
     "F_UBii = +/-F_rel*((L_GW_merg)/E_LEP)*Q_wave;"
     " L_GW = (32/5)*(G^4/c^5)*(m1*m2*(m1+m2)/r^5);"
     " tau_merge = 12/85*(a^4/G^3*m1*m2*(m1+m2))",
     "GW merger power and coalescence time (quadrupole formula full form)",
     "Peters 1964; a=semi-major axis; tau shrinks as a^4; GW inspiral",
     0.0, "W", "GW merger luminosity"),
]


def _esc(s: str) -> str:
    """Escape double-quotes for C++ string literals."""
    return s.replace("\\", "\\\\").replace('"', '\\"')


def get_section_E() -> str:
    lines = [
        "\n// ======================================================\n",
        "// SECTION E:  79 unique F_UBii variant types\n",
        "// Source: BB_C_Equations_04Sept2025.pdf\n",
        "//         (extracted from grok_share_c020496d9e.txt ~lines 692-6955)\n",
        "// General form:\n",
        "//   F_UBii,TYPE = +/-F_rel * (physics_factor/E_LEP) * Q_wave * scaling\n",
        "//   F_rel = 4.31e33 N  (2024 LEP)  Q_wave_mean = 3.97e4 J/m^3\n",
        "// ======================================================\n",
        "std::vector<BuoyancyEquation> buildSectionE() {\n",
        "    return {\n",
    ]
    for entry in _E_DATA:
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
            + '"' + _esc(system) + '", "E" },\n'
        )
    lines.append("    };\n}\n")
    return "".join(lines)
