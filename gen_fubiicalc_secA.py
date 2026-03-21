"""
gen_fubiicalc_secA.py
Returns C++ buildSectionA() — 29 per-system g_UQFF equations.

UQFF backbone shared by ALL 29 systems:
  g_UQFF(r,t) = (G*M(t)/r^2)*(1+H(t,z))*(1-B(t)/B_crit)*(1+F_env(t))
              + (Ug1+Ug2+Ug3+Ug4) + Lambda*c^2/3
              + (hbar/sqrt(dx*dp))*Int[psi*H*psi dV]*(2*pi/t_Hub)
              + rho_fluid*V*g + (M_vis+M_DM)*(d_rho/rho + 3GM/r^3)
              + SYSTEM_SPECIFIC_TERM
  H(t,z) = H0*sqrt(0.3*(1+z)^3 + 0.7)
"""

# Backbone shared by all systems (abbreviated for readability in expressions)
_BB = ("(G*M(t)/r^2)*(1+H(t,z))*(1-B/B_crit)*(1+F_env)"
       "+Ug1+Ug2+Ug3+Ug4+Lambda*c^2/3"
       "+(hbar/sqrt(dx*dp))*Int[psi*H*psi]*(2pi/t_Hub)"
       "+rho_f*V*g+(M_vis+M_DM)*(drho/rho+3GM/r^3)")

# (id, name, unique_term, context, note, system)
_SECTION_A_DATA = [
    (1,  "g_student",
     _BB + " [BASE: G*M_sun(t)/r(t)^2*(1+H0*t)]",
     "Student Guide Universe - cosmological base reference system",
     "Doc 1: calibration anchor; H0=67.4 km/s/Mpc; D_univ~92.77 Gly",
     "Student Guide Universe"),

    (2,  "g_magnetar_SGR1745",
     _BB + " + M_mag + D(t) + G*M_BH/r_BH^2",
     "Magnetar SGR 1745-2900 near Sgr A* - extreme B-field ~10^15 G",
     "Doc 2: M_mag=magnetic dipole torque; D(t)=drag; B_crit=4.4e13 T",
     "SGR 1745-2900 Magnetar"),

    (3,  "g_sgrA_star",
     _BB + " * sin(30deg) + G*M(t)^2/(c^4*r)*(dOmega/dt)^2",
     "Sagittarius A* SMBH - 4e6 Msun; galactic center accretion",
     "Doc 3: Kerr frame-dragging; a_spin~0.9; GWTC-4 compatible",
     "Sagittarius A*"),

    (4,  "g_tapestry",
     _BB + " + rho*v_wind^2",
     "Tapestry Starbirth Region - molecular cloud star formation; stellar winds",
     "Doc 4: v_wind~2000 km/s; SFR~10 Msun/yr",
     "Tapestry Starbirth Region"),

    (5,  "g_westerlund1",
     _BB + " + rho*v_wind^2 [cluster winds; implied]",
     "Westerlund 1 - massive young star cluster; ~10^5 Msun within 1 pc",
     "Doc 5: dense cluster; wind-feedback dominated",
     "Westerlund 1"),

    (6,  "g_westerlund2",
     _BB + " + rho*v_wind^2",
     "Westerlund 2 - star-forming giant HII region; eta Carinae neighbourhood",
     "Doc 6: FU_g1~2.43e-40 N; R(t)~-2.29e-41 N; FU_Bi~6.14e-32 N",
     "Westerlund 2"),

    (7,  "g_pillars",
     _BB + " * (1-E(t)) + rho*v_wind^2",
     "Pillars of Creation (M16 Eagle Nebula) - photoevaporation pillars",
     "Doc 7: E(t)=evaporation fraction; FU_g1~3.95e-41 N; FU_Bi~9.79e-33 N",
     "Pillars of Creation"),

    (8,  "g_rings_relativity",
     _BB + " * (1+L(t))",
     "Rings of Relativity - gravitational lensing system; Einstein ring",
     "Doc 8: L(t)=lensing amplification factor; arc geometry",
     "Rings of Relativity"),

    (9,  "g_ngc2525",
     _BB + " + G*M_BH/r_BH^2 - M_SN(t)",
     "NGC 2525 - SN Ia host galaxy; tidal interaction; AGN feedback",
     "Doc 10: M_SN(t)=supernova mass-loss term; Hubble SN Ia 2018",
     "NGC 2525"),

    (10, "g_ngc3603",
     _BB + " * (1-P(t)) + rho*v_wind^2",
     "NGC 3603 - most luminous HII region in Milky Way; O/WR star cluster",
     "Doc 11: P(t)=pressure ionization fraction",
     "NGC 3603"),

    (11, "g_bubble_nebula",
     _BB + " * (1+E(t)) + rho*v_exp^2",
     "Bubble Nebula (NGC 7635) - stellar wind blown bubble around O-star",
     "Doc 12: E(t)=expansion enhancement; v_exp~40 km/s",
     "Bubble Nebula NGC 7635"),

    (12, "g_antennae",
     _BB + " * (1-M_coll(t)) + rho*v_sf^2",
     "Antennae Galaxies (NGC 4038/4039) - ongoing major merger; starburst",
     "Doc 14: M_coll(t)=collision deceleration; SFR~20 Msun/yr",
     "Antennae Galaxies"),

    (13, "g_horsehead",
     _BB + " * (1-E(t)) + P_rad",
     "Horsehead Nebula (B33) - dark nebula silhouette against IC 434",
     "Doc 15: P_rad=radiation pressure from sigma Orionis",
     "Horsehead Nebula B33"),

    (14, "g_ngc1275",
     _BB + " + F_BH + M_fil",
     "NGC 1275 Perseus Cluster - brightest cluster galaxy; AGN jet filaments",
     "Doc 16: F_BH=AGN feedback force; M_fil=filament infall momentum",
     "NGC 1275 Perseus Cluster"),

    (15, "g_hubble_udf",
     _BB + " * (1+M_evo(t)) * (1-M_merge(t))",
     "Hubble Ultra Deep Field - z=0-7 galaxy population; cosmic web",
     "Doc 18: M_evo(t)=evolution modifier; M_merge(t)=merger fraction",
     "Hubble Ultra Deep Field"),

    (16, "g_ngc1792",
     _BB + " * (1+M_sf(t)) + F_sn",
     "NGC 1792 - active star-forming spiral; supernovae feedback",
     "Doc 19: M_sf(t)=SF modulator; F_sn=supernova kick force",
     "NGC 1792"),

    (17, "g_sombrero",
     _BB + " + G*M_BH/r_BH^2 + D_dust",
     "Sombrero Galaxy (M104) - prominent dust lane; 10^9 Msun SMBH",
     "Doc 20: D_dust=dust-lane pressure; bulge-disk decomposition",
     "Sombrero Galaxy M104"),

    (18, "g_saturn",
     "G*M_Sun/r_orbit^2 + T_ring + F_wind + " + _BB,
     "Saturn system - ring dynamics; magnetosphere; Enceladus plumes",
     "Doc 22: T_ring=ring torque; F_wind=solar wind; Cassini data",
     "Saturn System"),

    (19, "g_m16_eagle",
     _BB + " * (1+M_sf(t)) - E_rad",
     "M16 Eagle Nebula - young open cluster; Pillars host nebula",
     "Doc 23: E_rad=EUV erosion rate from cluster stars",
     "M16 Eagle Nebula"),

    (20, "g_crab_nebula",
     "G*M_NS/r(t)^2 + F_wind + M_mag + " + _BB,
     "Crab Nebula (M1) - pulsar wind nebula; 33 ms pulsar; 1054 CE SN",
     "Doc 24: M_mag=pulsar magnetic torque; F_wind=pulsar nebula wind; dE/dt~5e38 erg/s",
     "Crab Nebula M1"),

    (21, "g_hydrogen_atom",
     _BB + " * (1+P_term) + F_tech",
     "Hydrogen Atom - quantum scale; proton-electron Coulomb + nuclear",
     "Doc 27: P_term=pairing energy; F_tech=technology coupling term",
     "Hydrogen Atom"),

    (22, "g_hydrogen_resonance",
     "H_res=A_res*sin(2pi*f_res*t)+U_dp*SC_m*k_nuc+S_shell; "
     "A_res=k_A*Z*(A/A_H)*(1+d_pair); f_res=(E_bind/h)*(A_H/A)*(1+S_shell); "
     "U_dp=k*(A1*A2/f_dp^2)*cos(phi_dp); k_nuc=k0*(N/Z)*(1+d_pair); "
     "S_shell=0.1*(Z_magic+N_magic)",
     "Hydrogen Resonance - 7 sub-equations; nuclear + shell resonance",
     "Doc 28: SC_m~1; 7 coupled sub-equations for resonance frequency",
     "Hydrogen Resonance"),

    (23, "g_universe_diameter",
     "D_univ=2*D_p*(1+H(z)*t0)*(1+Lambda*c^2/(3*H0^2))"
     "*(1+(hbar/sqrt(dx*dp))*Int[psi*H*psi]/(G*M_tot))*(1+k*r_c^2)",
     "Universe Diameter - full UQFF cosmological scale equation",
     "Doc 26/29: D_univ~92.77 Gly (matches ~93 Gly observed); H0=67.4 km/s/Mpc",
     "Observable Universe"),

    (24, "g_ngc2264",
     _BB + " * (1+M_sf(t)) + delta_SFR",
     "NGC 2264 - Christmas Tree Cluster + Cone Nebula; pre-main-sequence stars",
     "Doc 18b / SOURCE115: 26D polynomial node; SFR~0.5 Msun/yr",
     "NGC 2264"),

    (25, "g_tadpole",
     _BB + " + G*dM_tidal/r_tidal^2",
     "Tadpole Galaxy (UGC 10214) - 280 kpc tidal tail from z~0.3 companion",
     "SOURCE115: tidal disruption; 26D resonance node",
     "Tadpole Galaxy UGC 10214"),

    (26, "g_mice_galaxies",
     _BB + " * (1-M_coll(t)) + rho*v_tidal^2",
     "Mice Galaxies (NGC 4676 A/B) - equal-mass merger; counter-rotating tails",
     "SOURCE115: merger-induced starburst; 26D polynomial master",
     "Mice Galaxies NGC 4676"),

    (27, "g_m42_orion",
     _BB + " + rho_PDR*v_shock^2 + P_OB_rad",
     "M42 Orion Nebula - photodissociation region; Trapezium OB cluster",
     "SOURCE115: PDR driven; Chandra X-ray data; 26D node",
     "M42 Orion Nebula"),

    (28, "g_carina_nebula",
     _BB + " + rho*v_wind^2 + F_eta_car",
     "Carina Nebula - eta Carinae LBV; multiple HII regions; ~60 O-stars",
     "SOURCE115: F_eta_car=LBV Great Eruption analogue; SFR~3 Msun/yr",
     "Carina Nebula"),

    (29, "g_eso137",
     _BB + " - rho_ICM*v_ram^2",
     "ESO 137-001 - ram pressure stripping in Norma cluster; x-ray tail",
     "SOURCE115: v_ram~2000 km/s; ICM rho~1e-27 kg/m^3; Chandra 2025",
     "ESO 137-001"),
]


def _esc(s: str) -> str:
    """Escape double-quotes for use inside C++ string literals."""
    return s.replace('"', '\\"')


def get_section_A() -> str:
    lines = [
        "\n// ======================================================\n",
        "// SECTION A:  29 per-system g_UQFF equations\n",
        "// UQFF backbone (ALL systems):\n",
        "//   g_UQFF(r,t) = (G*M(t)/r^2)*(1+H(t,z))*(1-B/B_crit)*(1+F_env)\n",
        "//               + (Ug1+Ug2+Ug3+Ug4) + Lambda*c^2/3\n",
        "//               + (hbar/sqrt(dx*dp))*Int[psi*H*psi dV]*(2pi/t_Hub)\n",
        "//               + rho_f*V*g + (M_vis+M_DM)*(drho/rho + 3GM/r^3)\n",
        "//               + SYSTEM_SPECIFIC_TERM\n",
        "//   H(t,z) = H0*sqrt(0.3*(1+z)^3 + 0.7)\n",
        "// ======================================================\n",
        "std::vector<BuoyancyEquation> buildSectionA() {\n",
        "    return {\n",
    ]
    for entry in _SECTION_A_DATA:
        id_, name, expr, ctx, note, system = entry
        lines.append(
            "        { " + str(id_) + ", "
            + '"' + _esc(name) + '", '
            + '"' + _esc(expr) + '", '
            + '"' + _esc(ctx)  + '", '
            + '"' + _esc(note) + '", '
            + "0.0, "
            + '"m/s^2", '
            + '"' + _esc(system) + '", "A" },\n'
        )
    lines.append("    };\n}\n")
    return "".join(lines)
