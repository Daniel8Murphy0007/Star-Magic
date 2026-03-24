// observational_systems_config.h
// System parameter configurations for observational astrophysical systems
// Use with existing SOURCE10 buoyancy physics terms in MAIN_1_CoAnQi.cpp
// Copyright - Daniel T. Murphy, 2025

#ifndef OBSERVATIONAL_SYSTEMS_CONFIG_H
#define OBSERVATIONAL_SYSTEMS_CONFIG_H

#include <map>
#include <string>

struct ObservationalSystem
{
    std::string name;
    std::string description;

    // Physical parameters
    double M;       // Mass (kg)
    double r;       // Radius (m)
    double L_X;     // X-ray luminosity (W)
    double B0;      // Magnetic field (T)
    double rho_gas; // Gas density (kg/m³)
    double T_gas;   // Gas temperature (K)
    double omega0;  // Angular frequency (rad/s)
    double t_age;   // System age/timescale (s)

    // Observational context
    std::string category;  // "galaxy_cluster", "pulsar", "agn", "nebula", "tde", "snr"
    std::string telescope; // Primary observation source
    
    // Epoch Framework Context (Grok Thread 4e0ecf23 - March 4, 2026)
    // Which cosmic epochs this system represents in Inflation/Force Chart
    // Epochs: 1=Fisile, 2=Star/Planet, 3=Galaxy/Quasar, 4=Magnetar/SMBH, 5=Globular
    int dominant_epoch;  // Primary epoch (1-5)
    bool epoch_1_present;  // Fisile Nuclei/Nebular (pre-stellar)
    bool epoch_2_present;  // Star/Planetary Atom (Ug1-3 active)
    bool epoch_3_present;  // Galaxies/Quasar (early Ug4)
    bool epoch_4_present;  // Magnetar/SMBH (Ug4 dominance)
    bool epoch_5_present;  // Globular Clusters (stabilization)
};

// ============================================================================
// SYSTEM DEFINITIONS (35+ observational targets from Source155-161)
// ============================================================================

static const std::map<std::string, ObservationalSystem> OBSERVATIONAL_SYSTEMS = {

    // ========== GALAXY CLUSTERS (Source155, 156, 157, 158, 159, 160) ==========

    {"ESO137", {"ESO 137-001", "Ram pressure stripping galaxy in Abell 3627 cluster",
                2e41,    // M (kg)
                6.17e21, // r (m)
                1e34,    // L_X (W)
                2e-9,    // B0 (T)
                1e-23,   // rho_gas (kg/m³)
                1e7,     // T_gas (K)
                1e-15,   // omega0 (rad/s)
                7.72e14, // t_age (s, ~24.5 Myr)
                "galaxy_cluster", "Chandra/ALMA/MUSE",
                3,       // Dominant Epoch: 3 (Galaxies/Quasar)
                false,   // No Epoch 1 (pre-stellar)
                true,    // Epoch 2 (stars in galaxy)
                true,    // Epoch 3 (galaxy structure) ← DOMINANT
                false,   // No Epoch 4 (no SMBH/magnetar signature)
                false}}, // No Epoch 5 (not globular cluster)

    {"NGC1365", {"NGC 1365", "Barred spiral galaxy with active nucleus",
                 7.17e41, // M (kg)
                 9.46e20, // r (m)
                 1e36,    // L_X (W)
                 1e-9,    // B0 (T)
                 1e-23,   // rho_gas (kg/m³)
                 2e7,     // T_gas (K)
                 1e-15,   // omega0 (rad/s)
                 1.1e16,  // t_age (s, ~350 Myr)
                 "agn", "Chandra/JWST"}},

    {"ElGordo", {"El Gordo (ACT-CL J0102-4915)", "Most massive distant galaxy cluster merger",
                 4.97e45, // M (kg)
                 3.09e22, // r (m)
                 2e38,    // L_X (W)
                 1e-10,   // B0 (T)
                 1e-24,   // rho_gas (kg/m³)
                 1.4e8,   // T_gas (K)
                 1e-15,   // omega0 (rad/s)
                 2.21e16, // t_age (s, ~700 Myr)
                 "galaxy_cluster", "Chandra/ALMA/VLT"}},

    {"PLCKG287", {"PLCK G287.0+32.9", "Planck-detected galaxy cluster",
                  1.989e44, // M (kg)
                  3.09e22,  // r (m)
                  5e37,     // L_X (W)
                  5e-10,    // B0 (T)
                  5e-24,    // rho_gas (kg/m³)
                  8e7,      // T_gas (K)
                  1e-15,    // omega0 (rad/s)
                  1.5e16,   // t_age (s)
                  "galaxy_cluster", "Planck/Chandra"}},

    {"PSZ2G181", {"PSZ2 G181.06+48.47", "Planck SZ2 catalog galaxy cluster",
                  1.989e44, // M (kg)
                  3.09e22,  // r (m)
                  4e37,     // L_X (W)
                  4e-10,    // B0 (T)
                  4e-24,    // rho_gas (kg/m³)
                  7e7,      // T_gas (K)
                  1e-15,    // omega0 (rad/s)
                  1.8e16,   // t_age (s)
                  "galaxy_cluster", "Planck/XMM"}},

    {"NGC4839", {"NGC 4839", "Galaxy in Coma Cluster showing ram pressure effects",
                 3e41,   // M (kg)
                 8e20,   // r (m)
                 8e35,   // L_X (W)
                 1.5e-9, // B0 (T)
                 2e-23,  // rho_gas (kg/m³)
                 1.5e7,  // T_gas (K)
                 1e-15,  // omega0 (rad/s)
                 1e16,   // t_age (s)
                 "galaxy_cluster", "Chandra/HST"}},

    {"Abell2256", {"Abell 2256", "Merging galaxy cluster with radio relics",
                   1.23e45, // M (kg)
                   3.93e22, // r (m)
                   1.5e38,  // L_X (W)
                   2e-10,   // B0 (T)
                   1e-24,   // rho_gas (kg/m³)
                   9e7,     // T_gas (K)
                   1e-15,   // omega0 (rad/s)
                   2e16,    // t_age (s)
                   "galaxy_cluster", "Chandra/VLA"}},

    {"M84", {"M84 (NGC 4374)", "Lenticular galaxy in Virgo Cluster with jets",
             1.46e45, // M (kg)
             3.09e22, // r (m)
             3e37,    // L_X (W)
             8e-10,   // B0 (T)
             3e-24,   // rho_gas (kg/m³)
             6e7,     // T_gas (K)
             1e-15,   // omega0 (rad/s)
             1.2e16,  // t_age (s)
             "agn", "Chandra/HST/VLA"}},

    // ========== PULSARS (Source155, 156, 161) ==========

    {"Vela", {"Vela Pulsar", "Young pulsar with pulsar wind nebula",
              2.8e30,  // M (kg, ~1.4 solar masses)
              1.7e17,  // r (m, ~170 Gm standoff distance)
              1e27,    // L_X (W)
              3e-8,    // B0 (T, surface field ~1e8 G)
              1e-23,   // rho_gas (kg/m³)
              1e6,     // T_gas (K, PWN)
              1e-12,   // omega0 (rad/s, ~89 ms period)
              3.47e11, // t_age (s, ~11,000 years)
              "pulsar", "Chandra/Fermi",
              4,       // Dominant Epoch: 4 (Magnetar/SMBH scale)
              false,   // No Epoch 1
              true,    // Epoch 2 (formed from star)
              true,    // Epoch 3 (in galaxy)
              true,    // Epoch 4 (extreme field) ← DOMINANT
              false}}, // No Epoch 5

    {"J1610", {"J1610+1811", "Millisecond pulsar in globular cluster",
               2.785e30, // M (kg)
               3.09e15,  // r (m, light cylinder radius)
               5e26,     // L_X (W)
               1e-8,     // B0 (T)
               1e-22,    // rho_gas (kg/m³)
               5e5,      // T_gas (K)
               1e-11,    // omega0 (rad/s)
               1e15,     // t_age (s, old recycled)
               "pulsar", "Chandra/XMM"}},

    {"Crab", {"Crab Nebula Pulsar", "Young energetic pulsar with wind nebula",
              1e31,    // M (kg, includes nebula mass)
              4.73e16, // r (m, nebula radius)
              1e31,    // L_X (W, total nebula)
              5e-8,    // B0 (T, pulsar surface ~1e12 G)
              1e-22,   // rho_gas (kg/m³)
              1e7,     // T_gas (K)
              2e-10,   // omega0 (rad/s, 33 ms)
              2.9e10,  // t_age (s, ~900 years)
              "pulsar", "Chandra/HST/JWST"}},

    {"ASKAPJ1832", {"ASKAP J1832-0911", "Radio pulsar discovered by ASKAP",
                    2.785e30, // M (kg)
                    4.63e16,  // r (m)
                    3e26,     // L_X (W)
                    2e-8,     // B0 (T)
                    1e-23,    // rho_gas (kg/m³)
                    8e5,      // T_gas (K)
                    1e-12,    // omega0 (rad/s)
                    5e14,     // t_age (s)
                    "pulsar", "ASKAP/Parkes"}},

    // ========== ACTIVE GALACTIC NUCLEI (Source155, 157, 158) ==========

    {"CentaurusA", {"Centaurus A (NGC 5128)", "Nearest radio galaxy with prominent jets",
                    1.094e38, // M (kg, SMBH mass ~5.5e7 solar)
                    6.17e17,  // r (m, jet interaction region)
                    2e34,     // L_X (W)
                    1e-6,     // B0 (T, jet magnetic field)
                    1e-21,    // rho_gas (kg/m³)
                    5e6,      // T_gas (K)
                    1e-12,    // omega0 (rad/s)
                    1e14,     // t_age (s)
                    "agn", "Chandra/ALMA/VLT",
                    4,        // Dominant Epoch: 4 (SMBH)
                    false,    // No Epoch 1
                    true,     // Epoch 2 (stars in galaxy)
                    true,     // Epoch 3 (galaxy formed)
                    true,     // Epoch 4 (SMBH jets) ← DOMINANT
                    false}},  // No Epoch 5

    {"M104", {"M104 (Sombrero Galaxy)", "Edge-on spiral with large bulge and SMBH",
              1.5e38, // M (kg, SMBH ~8e8 solar)
              5e20,   // r (m)
              5e33,   // L_X (W)
              5e-10,  // B0 (T)
              5e-24,  // rho_gas (kg/m³)
              3e6,    // T_gas (K)
              1e-13,  // omega0 (rad/s)
              2e16,   // t_age (s)
              "agn", "Chandra/HST"}},

    // ========== TIDAL DISRUPTION EVENT (Source155) ==========

    {"ASASSN14li", {"ASASSN-14li", "Tidal disruption event with jetted outflow",
                    1.989e37, // M (kg, SMBH ~3e6 solar)
                    3.09e18,  // r (m, tidal radius)
                    1e37,     // L_X (W, peak)
                    1e-5,     // B0 (T, accretion disk field)
                    1e-21,    // rho_gas (kg/m³, debris)
                    1e7,      // T_gas (K)
                    1e-12,    // omega0 (rad/s)
                    9.504e6,  // t_age (s, ~110 days since disruption)
                    "tde", "Chandra/Swift/ASAS-SN"}},

    // ========== NEBULAE (Source157, 158, 159, 160) ==========

    {"NGC346", {"NGC 346", "Star-forming region in Small Magellanic Cloud",
                1e36,   // M (kg)
                1.5e17, // r (m, ~50 pc)
                2e32,   // L_X (W)
                1e-9,   // B0 (T)
                1e-20,  // rho_gas (kg/m³)
                1e4,    // T_gas (K, HII region)
                1e-14,  // omega0 (rad/s)
                1e14,   // t_age (s, ~3 Myr)
                "nebula", "HST/JWST/Chandra",
                2,      // Dominant Epoch: 2 (Star formation)
                false,  // No Epoch 1 (already past fisile)
                true,   // Epoch 2 (active star formation) ← DOMINANT
                false,  // No Epoch 3 (not galaxy-scale)
                false,  // No Epoch 4
                false}},// No Epoch 5

    {"M16", {"M16 (Eagle Nebula)", "Star-forming region with Pillars of Creation",
             1e36,    // M (kg)
             2.36e17, // r (m)
             5e32,    // L_X (W)
             2e-9,    // B0 (T)
             1e-20,   // rho_gas (kg/m³)
             8e3,     // T_gas (K)
             1e-14,   // omega0 (rad/s)
             2e14,    // t_age (s, ~5.5 Myr)
             "nebula", "HST/Spitzer/Chandra"}},

    {"NGC1672", {"NGC 1672", "Barred spiral with intense star formation",
                 5e41,   // M (kg)
                 7e20,   // r (m)
                 2e36,   // L_X (W)
                 1.2e-9, // B0 (T)
                 1e-23,  // rho_gas (kg/m³)
                 1.5e7,  // T_gas (K)
                 1e-15,  // omega0 (rad/s)
                 1.5e16, // t_age (s)
                 "nebula", "HST/Chandra"}},

    {"Tarantula", {"Tarantula Nebula (30 Doradus)", "Giant HII region in Large Magellanic Cloud",
                   1e36,  // M (kg)
                   2e17,  // r (m, ~65 pc)
                   1e33,  // L_X (W)
                   3e-9,  // B0 (T)
                   1e-20, // rho_gas (kg/m³)
                   1e4,   // T_gas (K)
                   1e-14, // omega0 (rad/s)
                   3e14,  // t_age (s, ~8 Myr)
                   "nebula", "HST/Spitzer/Chandra"}},

    // ========== SUPERNOVA REMNANTS (Source158, 159, 160) ==========

    {"Tycho", {"Tycho's Supernova Remnant", "Type Ia supernova remnant from 1572",
               1e31,    // M (kg, ejecta + swept-up)
               1e17,    // r (m, ~3 pc)
               5e30,    // L_X (W)
               1e-7,    // B0 (T, shock magnetic field)
               1e-22,   // rho_gas (kg/m³)
               2e7,     // T_gas (K, shock-heated)
               1e-13,   // omega0 (rad/s)
               1.43e10, // t_age (s, ~453 years)
               "snr", "Chandra/XMM/VLA"}},

    // SNR G272.2-03.2 — Added Session 133 (grok_share_c35c3b7a1, Nov 24-28 2025)
    // Chandra "cosmic gourd" Nov 2025 release; Type Ia thermal composite; Vela constellation
    // Cosmic Egg connection: egg-hatching shell morphology; SCm shell collapse (Type Ia no NS)
    {"SNR_G272", {"SNR G272.2-03.2 (Cosmic Gourd)", "Type Ia supernova remnant in Vela; thermal composite; Chandra Nov 2025 release; ~7500 yr; UQFF egg-hatching structure",
                  4e30,    // M (kg, ~2 M_sun ejecta + swept ISM at 7500 yr)
                  4e17,    // r (m, ~13 pc at 2.5 kpc, 17.5 arcmin diameter)
                  2e29,    // L_X (W, cooling older remnant, Chandra detected)
                  5e-8,    // B0 (T, compressed ISM magnetic field at shock)
                  5e-23,   // rho_gas (kg/m³, Vela ISM density)
                  1e7,     // T_gas (K, shock-heated, cooling at 7500 yr)
                  0.0,     // omega0 (rad/s, Type Ia — no neutron star remnant)
                  2.37e11, // t_age (s, ~7500 years × 3.156e7 s/yr)
                  "snr", "Chandra/XMM-Newton",
                  2,       // Dominant Epoch: 2 (Star/Planetary — SN remnant post-stellar death)
                  false,   // No Epoch 1 (pre-stellar phase past)
                  true,    // Epoch 2 (stellar death/remnant expansion) ← DOMINANT
                  false,   // No Epoch 3 (not galactic scale)
                  false,   // No Epoch 4 (no SMBH/magnetar — Type Ia)
                  false}}, // No Epoch 5 (not globular cluster)

    // ========== GALAXIES (Source157, 158, 160) ==========

    {"M74", {"M74 (NGC 628)", "Grand design spiral galaxy",
             7.17e41, // M (kg)
             9.46e20, // r (m)
             5e35,    // L_X (W)
             8e-10,   // B0 (T)
             1e-23,   // rho_gas (kg/m³)
             1e7,     // T_gas (K)
             1e-15,   // omega0 (rad/s)
             2e16,    // t_age (s)
             "galaxy", "HST/Chandra/JWST"}},

    {"NGC253", {"NGC 253 (Sculptor Galaxy)", "Starburst galaxy with superwind",
                4e40,  // M (kg)
                4e20,  // r (m)
                1e36,  // L_X (W)
                2e-9,  // B0 (T)
                2e-23, // rho_gas (kg/m³)
                2e7,   // T_gas (K, hot wind)
                1e-15, // omega0 (rad/s)
                1e16,  // t_age (s)
                "galaxy", "Chandra/ALMA/VLT"}},

    // ========== MULTI-WAVELENGTH COLLECTIONS (Source156, 161) ==========

    {"Sonification", {"Chandra Sonification Collection", "Multi-system sonified X-ray data collection",
                      1.989e31, // M (kg, representative)
                      6.17e16,  // r (m, representative)
                      1e32,     // L_X (W, representative)
                      1e-9,     // B0 (T, representative)
                      1e-22,    // rho_gas (kg/m³)
                      1e6,      // T_gas (K)
                      1e-13,    // omega0 (rad/s)
                      1e15,     // t_age (s)
                      "multi_system", "Chandra"}},

    {"ChandraWebb", {"Chandra-Webb Collaborative Observations", "Combined X-ray and infrared multi-system survey",
                     1e40,  // M (kg, representative)
                     1e20,  // r (m, representative)
                     1e35,  // L_X (W, representative)
                     1e-9,  // B0 (T, representative)
                     1e-23, // rho_gas (kg/m³)
                     1e7,   // T_gas (K)
                     1e-14, // omega0 (rad/s)
                     1e16,  // t_age (s)
                     "multi_system", "Chandra/JWST"}},

    {"SupernovaSurvey", {"Supernova Observational Survey", "Multi-supernova comparative study",
                         1e30,  // M (kg, representative NS)
                         1e10,  // r (m, representative remnant)
                         1e31,  // L_X (W, representative)
                         1e-7,  // B0 (T, representative)
                         1e-22, // rho_gas (kg/m³)
                         1e7,   // T_gas (K)
                         1e-13, // omega0 (rad/s)
                         1e11,  // t_age (s, representative ~3000 yrs)
                         "multi_system", "Chandra/HST/Swift"}},

    // =========================================================================
    // GROK THREAD 0904a12a — 5 NEW SYSTEMS (March 6, 2026)
    // Source: GrokThread_UQFF_0904_Validation.py (systems 25-29)
    // =========================================================================

    {"GRO_J1655-40", {"GRO J1655-40 Micro-quasar", "Black hole X-ray binary with superluminal jets; UQFF micro-quasar test case",
                      1.28e31,  // M (kg, ~6.4 M_sun BH)
                      1.5e9,    // r (m, ~2.2 R_sun companion)
                      1e30,     // L_X (W, near-Eddington outburst)
                      1e-5,     // B0 (T, jet field estimate)
                      1e-19,    // rho_gas (kg/m³, accretion disc)
                      1e7,      // T_gas (K, disc corona)
                      4.4e-5,   // omega0 (rad/s, 2.62-day orbital period)
                      3.15e14,  // t_age (s, ~10 Myr)
                      "micro_quasar", "RXTE/ASCA/HST"}},

    {"CygnusLoop", {"Cygnus Loop (Veil Nebula)", "Middle-aged SNR ~20,000 yr; shock-heated ISM; UQFF SNR test case",
                    2.78e30,  // M (kg, ~1.4 M_sun ejecta + swept ISM)
                    1.9e19,   // r (m, ~2 Mpc shock radius)
                    1e30,     // L_X (W, X-ray luminosity)
                    1e-8,     // B0 (T, post-shock field)
                    5e-24,    // rho_gas (kg/m³, ISM density)
                    3e6,      // T_gas (K, shock temperature)
                    0.0,      // omega0 (rad/s, no rotation)
                    6.3e11,   // t_age (s, ~20,000 yr)
                    "snr", "ROSAT/XMM/Chandra"}},

    {"G292.0+1.8", {"G292.0+1.8 SNR/PWN", "Young oxygen-rich SNR with central pulsar; UQFF PWN/SNR test case",
                    2.78e30,  // M (kg, ~1.4 M_sun NS + ejecta)
                    1.5e17,   // r (m, ~5 pc radius)
                    3e30,     // L_X (W, X-ray luminosity)
                    1e-9,     // B0 (T, pulsar wind nebula field)
                    1e-22,    // rho_gas (kg/m³, ejecta density)
                    1e7,      // T_gas (K, shock temperature)
                    2.9e1,    // omega0 (rad/s, ~135 ms pulsar)
                    4.7e10,   // t_age (s, ~1,500 yr)
                    "snr_pwn", "Chandra/XMM"}},

    {"NGC7293", {"NGC 7293 Helix Nebula", "Nearest planetary nebula; UQFF PN vacuum density test",
                 1.19e30,  // M (kg, ~0.6 M_sun white dwarf)
                 1.54e15,  // r (m, ~0.5 pc inner shell)
                 1e26,     // L_X (W, UV/X-ray from hot WD)
                 1e-7,     // B0 (T, WD surface field estimate)
                 1e-21,    // rho_gas (kg/m³, nebular density)
                 1e4,      // T_gas (K, nebular temperature)
                 0.0,      // omega0 (rad/s)
                 6.3e12,   // t_age (s, ~200 kyr)
                 "planetary_nebula", "Chandra/HST/Spitzer"}},

    {"PerseusCluster", {"Perseus Galaxy Cluster (Abell 426)", "Brightest X-ray cluster; giant cavities from AGN jets; UQFF ICM test",
                        6.65e44,  // M (kg, ~3.3e14 M_sun total)
                        2.31e22,  // r (m, ~750 kpc virial radius)
                        7e37,     // L_X (W, X-ray luminosity)
                        1e-9,     // B0 (T, ICM magnetic field)
                        1e-25,    // rho_gas (kg/m³, ICM core density)
                        7e7,      // T_gas (K, ICM temperature ~6 keV)
                        0.0,      // omega0 (rad/s)
                        4.1e17,   // t_age (s, ~13 Gyr)
                        "galaxy_cluster", "Chandra/XMM/Hitomi"}},

    // =========================================================================
    // PAPER_430-445: Per-System MUGE Library (Session 119, grok_share_68eb34022.txt)
    // Author: Daniel T. Murphy  |  Integration Date: 2026-03-22
    // =========================================================================

    // PAPER_430 — SGR 0501+4516 Magnetar (first complete per-system MUGE; B(t) decay)
    {"SGR0501_4516", {"SGR 0501+4516 Magnetar", "First complete per-system MUGE; B(t) exponential decay; anomalous separation from HB9 SNR (PAPER_430)",
                      2.785e30,  // M (kg, 1.4 M_sun)
                      2.0e4,     // r (m, 20 km radius)
                      1e31,      // L_X (W, X-ray luminosity)
                      1e10,      // B0 (T, magnetar surface field)
                      1e-20,     // rho_gas (kg/m³, magnetosphere)
                      1e8,       // T_gas (K, surface temperature)
                      1.222,     // omega0 (rad/s, rotation ~5.76 s period)
                      1.578e11,  // t_age (s, ~5000 yr)
                      "magnetar", "Chandra/XMM/NICER"}},

    // PAPER_431 — SGR 1745-2900 near Sgr A* (B(t) decay + BH proximity gravitational coupling)
    {"SGR1745_2900", {"SGR 1745-2900 Magnetar near Sgr A*", "Near-GC magnetar; B(t) decay + g_BH tidal coupling from Sgr A*; cumulative decay energy (PAPER_431)",
                      2.785e30,  // M (kg, 1.4 M_sun)
                      1.0e4,     // r (m, 10 km radius)
                      4e27,      // L_X (W, L0 initial luminosity)
                      2e10,      // B0 (T, stronger field than SGR0501)
                      1e-20,     // rho_gas (kg/m³, magnetosphere)
                      1e8,       // T_gas (K)
                      0.5,       // omega0 (rad/s)
                      3.156e10,  // t_age (s, ~1000 yr)
                      "magnetar", "Chandra/XMM/NuSTAR"}},

    // PAPER_432 — Sagittarius A* SMBH (M(t) accretion + DM precession sin(30°))
    {"SgrA_SMBH", {"Sagittarius A* SMBH", "GC SMBH; M(t) accretion model; dark matter precession term sin(30 deg); 9 Gyr evolution (PAPER_432)",
                   8.553e36,  // M (kg, 4.3e6 M_sun)
                   1.27e10,   // r (m, ~3 Schwarzschild radii)
                   1e37,      // L_X (W, Sgr A* quiescent luminosity)
                   1.0,       // B0 (T, near-horizon field estimate)
                   1e-19,     // rho_gas (kg/m³, accretion disk)
                   1e10,      // T_gas (K, accretion torus)
                   0.3,       // omega0 (dimensionless spin parameter a)
                   2.84e17,   // t_age (s, ~9 Gyr)
                   "black_hole", "Chandra/VLA/EHT/Gravity"}},

    // PAPER_433 — Tapestry of Blazing Starbirth LMC (M(t) SF wind feedback)
    {"TapestryStarbirth", {"Tapestry of Blazing Starbirth (LMC N11)", "LMC star-forming region; M(t) wind/SF feedback; multi-cluster OB associations (PAPER_433)",
                           4.774e32,  // M (kg, ~2.4e2 M_sun total stellar)
                           9.461e16,  // r (m, ~10 ly half-radius)
                           1e33,      // L_X (W, UV/optical luminosity)
                           1e-6,      // B0 (T, ISM magnetic field)
                           1e-21,     // rho_gas (kg/m³, HII region density)
                           1e4,       // T_gas (K, HII region ~10^4 K)
                           0.0,       // omega0 (rad/s)
                           1.578e14,  // t_age (s, ~5 Myr)
                           "star_forming", "HST/Spitzer/Chandra"}},

    // PAPER_434 — Westerlund 2 YMC (M(t) growth; τ_SF=2 Myr; v_wind=2000 km/s)
    {"Westerlund2", {"Westerlund 2 Young Massive Cluster", "Most massive YMC in MW; M0=30000 M_sun; Mf=3.333 growth; supersonic wind 2000 km/s; tau_SF=2 Myr (PAPER_434)",
                     5.967e34,  // M (kg, 30000 M_sun initial)
                     9.461e16,  // r (m, ~10 ly half-radius)
                     1e34,      // L_X (W, UV bolometric luminosity)
                     1e-5,      // B0 (T, wind-swept ISM)
                     1e-20,     // rho_gas (kg/m³, wind density)
                     3e4,       // T_gas (K, stellar wind shocked gas)
                     0.0,       // omega0 (rad/s)
                     6.31e13,   // t_age (s, ~2 Myr)
                     "star_cluster", "Chandra/HST/VLT"}},

    // PAPER_435 — Pillars of Creation M16 NGC 6611 (erosion coupling E(t); τ_erosion=1 Myr)
    {"PillarsOfCreation", {"Pillars of Creation (Eagle Nebula M16, NGC 6611)", "Iconic photo-eroding molecular pillars; E(t)=0.1 exp(-t/tau) gravity suppression; tau=1 Myr (PAPER_435)",
                           2.009e34,  // M (kg, 10100 M_sun initial)
                           4.731e16,  // r (m, ~5 ly pillar half-length)
                           1e33,      // L_X (W, UV ionization luminosity)
                           1e-6,      // B0 (T, pillar magnetic field)
                           1e-21,     // rho_gas (kg/m³, molecular cloud density)
                           1e4,       // T_gas (K, photo-ionized gas)
                           0.0,       // omega0 (rad/s)
                           3.156e13,  // t_age (s, ~1 Myr)
                           "nebula", "HST/JWST/Spitzer"}},

    // PAPER_436 — Rings of Relativity GAL-CLUS-022058s Einstein ring (z=0.5; lensing amplification L)
    {"RingsOfRelativity", {"Rings of Relativity (GAL-CLUS-022058s Molten Ring)", "Near-perfect Einstein ring at z=0.5; lensing amplification factor L=(GM/c^2 r_E)(D_LS/D_S) MUGE correction (PAPER_436)",
                           1.989e44,  // M (kg, 1e14 M_sun lensing cluster)
                           3.086e20,  // r (m, 10 kpc Einstein radius)
                           1e38,      // L_X (W, cluster X-ray luminosity)
                           1e-5,      // B0 (T, ICM magnetic field)
                           1e-21,     // rho_gas (kg/m³, ICM density)
                           5e7,       // T_gas (K, ICM temperature)
                           0.0,       // omega0 (rad/s)
                           1.262e17,  // t_age (s, ~4 Gyr at z=0.5)
                           "gravitational_lens", "Hubble/WFC3"}},

    // PAPER_438 — NGC 2525 Galaxy (SN 2018gv host; SN mass loss + g_BH proximity)
    {"NGC2525", {"Galaxy NGC 2525 (SN 2018gv host)", "SN host galaxy; supernova mass loss channel M_SN=8 M_sun; BH proximity g_BH=GMbh/r_BH^2; z=0.016 (PAPER_438)",
                 1.989e40,  // M (kg, ~1e10 M_sun galaxy mass)
                 2.836e20,  // r (m, ~9 kpc disk radius)
                 1e36,      // L_X (W, galaxy bolometric)
                 1e-5,      // B0 (T, ISM magnetic field)
                 1e-25,     // rho_gas (kg/m³, ISM density)
                 1e4,       // T_gas (K, ISM)
                 1e-16,     // omega0 (rad/s, galactic rotation)
                 4.73e17,   // t_age (s, ~15 Gyr)
                 "galaxy", "HST/Chandra"}},

    // PAPER_439 — NGC 3603 Extreme Star Cluster (cavity pressure P(t); dual stellar wind)
    {"NGC3603", {"NGC 3603 Extreme Star Cluster", "Most luminous compact star cluster in MW; P(t)=4e-8 exp(-t/tau) cavity pressure; dual O+WR wind channels (PAPER_439)",
                 7.956e35,  // M (kg, ~400000 M_sun total)
                 8.988e16,  // r (m, ~9.5 ly half-radius)
                 1e35,      // L_X (W, UV bolometric luminosity)
                 1e-5,      // B0 (T, wind-swept field)
                 1e-20,     // rho_gas (kg/m³, wind density)
                 3e4,       // T_gas (K)
                 0.0,       // omega0 (rad/s)
                 3.156e13,  // t_age (s, ~1 Myr)
                 "star_cluster", "HST/VLT/Chandra"}},

    // PAPER_440 — Bubble Nebula NGC 7635 (growing expansion E(t); OB wind bubble)
    {"BubbleNebula", {"Bubble Nebula NGC 7635", "OB-star wind-blown bubble; growing expansion term E(t)=E0*(1-exp(-t/tau)); 4 Myr age (PAPER_440)",
                      9.149e31,  // M (kg, ~46 M_sun central star)
                      4.731e16,  // r (m, ~5 ly bubble radius)
                      1e32,      // L_X (W, UV/wind luminosity)
                      1e-6,      // B0 (T, ISM magnetic field)
                      1e-21,     // rho_gas (kg/m³, swept-up shell density)
                      1e4,       // T_gas (K, ionized gas)
                      0.0,       // omega0 (rad/s)
                      1.262e14,  // t_age (s, ~4 Myr)
                      "nebula", "HST/Chandra/Spitzer"}},

    // PAPER_441 — Antennae Galaxies NGC 4038+4039 (merger interaction boost I(t))
    {"AntennaeGalaxies", {"Antennae Galaxies NGC 4038+4039", "Archetypal merging galaxy pair; merger boost I(t)=I0*exp(-t/tau_merge); SFR=20 M_sun/yr; z=0.0105 (PAPER_441)",
                          3.978e41,  // M (kg, ~2e11 M_sun combined)
                          2.838e20,  // r (m, ~9.2 kpc effective radius)
                          1e36,      // L_X (W, starburst X-ray)
                          1e-5,      // B0 (T, starburst ISM field)
                          1e-21,     // rho_gas (kg/m³, ISM)
                          1e4,       // T_gas (K, ISM)
                          5e-17,     // omega0 (rad/s, galactic rotation)
                          1.262e16,  // t_age (s, ~400 Myr merger)
                          "galaxy_merger", "HST/Chandra/VLA"}},

    // PAPER_442 — Horsehead Nebula Barnard 33 (growing erosion E(t); τ=5 Myr)
    {"HorseheadNebula", {"Horsehead Nebula (Barnard 33)", "Iconic photo-eroded molecular pillar; growing erosion E(t)=E0*(1-exp(-t/tau)); tau=5 Myr (PAPER_442)",
                         1.989e33,  // M (kg, ~1000 M_sun molecular cloud)
                         2.365e16,  // r (m, ~2.5 ly) 
                         1e31,      // L_X (W, illuminating star sigma Ori)
                         1e-6,      // B0 (T, molecular cloud field)
                         1e-21,     // rho_gas (kg/m³, molecular cloud)
                         50.0,      // T_gas (K, cold molecular gas)
                         0.0,       // omega0 (rad/s)
                         1.578e14,  // t_age (s, ~5 Myr erosion)
                         "nebula", "HST/JWST/ALMA"}},

    // PAPER_443 — NGC 1275 Perseus A BCG (B(t) decay + F(t) filament coupling + cooling flow)
    {"NGC1275", {"NGC 1275 Perseus A BCG", "Perseus cluster BCG; B(t) decay + F(t)=F0*exp(-t/tau) Halpha filament + ICM cooling flow channel (PAPER_443)",
                 1.989e41,  // M (kg, ~1e11 M_sun stellar)
                 1.892e21,  // r (m, ~61 kpc effective radius)
                 5e37,      // L_X (W, X-ray cooling flow luminosity)
                 5e-9,      // B0 (T, ICM magnetic field from Faraday rotation)
                 1e-20,     // rho_gas (kg/m³, ICM core density)
                 3e7,       // T_gas (K, ICM ~2.5 keV)
                 0.0,       // omega0 (rad/s)
                 3.156e15,  // t_age (s, ~100 Myr cooling flow age)
                 "bcg", "Chandra/XMM/Hitomi/HST"}},

    // PAPER_444 — Hubble Ultra Deep Field z=3.5 (cosmic scale high-redshift MUGE)
    {"HUDFGalaxies", {"Hubble Ultra Deep Field (z=3.5 ensemble)", "Cosmological field; highest-z MUGE derivation; cosmic scale gravity at 3.5 <= z <= 3.5; Omega_m H^3 corrections (PAPER_444)",
                      1.989e42,  // M (kg, ~1e12 M_sun field ensemble)
                      1.230e27,  // r (m, ~40 Mpc co-moving scale)
                      1e38,      // L_X (W, UV/optical ensemble)
                      1e-10,     // B0 (T, intergalactic magnetic field)
                      1e-22,     // rho_gas (kg/m³, IGM density at z=3.5)
                      1e4,       // T_gas (K, reionization-era IGM)
                      0.0,       // omega0 (rad/s)
                      3.156e16,  // t_age (s, ~1 Gyr at z=3.5)
                      "cosmic_field", "HST/JWST/VLT"}},

    // PAPER_445 — NGC 1792 Stellar Forge (starburst wind dominance channel)
    {"NGC1792", {"NGC 1792 Stellar Forge", "Starburst galaxy; wind dominance channel T9_wind >> T1_gravity; SFR ~5 M_sun/yr; z=0.0095 (PAPER_445)",
                 1.989e40,  // M (kg, ~1e10 M_sun)
                 7.569e20,  // r (m, ~24.5 kpc disk radius)
                 1e36,      // L_X (W, starburst luminosity)
                 1e-5,      // B0 (T, starburst ISM field)
                 1e-21,     // rho_gas (kg/m³, ISM)
                 1e4,       // T_gas (K, ISM)
                 1e-16,     // omega0 (rad/s, galactic rotation)
                 3.156e15,  // t_age (s, ~100 Myr starburst)
                 "galaxy", "HST/Chandra/Spitzer"}},

    // PAPER_447 — M51 Whirlpool Galaxy (tidal NGC5195, 2-arm spiral density wave)
    {"M51Whirlpool", {"M51 Whirlpool Galaxy (+ NGC 5195)", "Interacting spiral; 2-arm spiral density wave psi=A*exp(-r^2/2sigma^2)*exp(i*(m*theta-omega*t)); tidal merger tau=5e8 yr; BH reaction Ug4; Session 120 (PAPER_447)",
                      3.182e41,  // M (kg, ~1.6e11 M_sun)
                      7.277e20,  // r (m, ~23.58 kpc)
                      1e37,      // L_X (W, spiral galaxy)
                      1e-5,      // B0 (T, ISM field)
                      1e-21,     // rho_gas (kg/m^3)
                      1e4,       // T_gas (K)
                      1e-16,     // omega0 (rad/s, galactic rotation)
                      3.156e16,  // t_age (s, ~1 Gyr)
                      "galaxy", "HST/Chandra/Spitzer/JWST"}},

    // PAPER_448 — V838 Monocerotis (expanding light echo, peak L=2.3e38 W)
    {"V838Mon", {"V838 Monocerotis (light echo nova)", "Peculiar stellar eruption; expanding light echo R_echo(t)=c*(t_obs-D/c); Ug_echo=G*M/R_echo^2*exp(-R_echo/r0); L_peak=2.3e38 W; Session 120 (PAPER_448)",
                 1.591e31,  // M (kg, ~8 M_sun progenitor)
                 5.677e16,  // r (m, ~6 ly light echo radius)
                 2.3e38,    // L_X (W, peak eruption luminosity)
                 1e-4,      // B0 (T, stellar field)
                 1e-22,     // rho_gas (kg/m^3, circumstellar)
                 3500.0,    // T_gas (K, cool giant)
                 1e-6,      // omega0 (rad/s, stellar rotation)
                 6.31e8,    // t_age (s, ~20 yr post-eruption)
                 "stellar", "HST/Spitzer/VLT/Keck"}},

    // PAPER_449 — NGC1300 Barred Galaxy (bar arm resonance v_arm=200 km/s)
    {"NGC1300Bar", {"NGC 1300 Barred Spiral Galaxy", "Barred spiral; Ug3_arm=G*M/r_arm^2 bar resonance; v_arm=200 km/s; orbital resonance forcing; Session 120 (PAPER_449)",
                    1.989e41,  // M (kg, ~1e11 M_sun)
                    3.638e20,  // r (m, ~11.79 kpc)
                    1e36,      // L_X (W)
                    1e-5,      // B0 (T)
                    1e-21,     // rho_gas (kg/m^3)
                    1e4,       // T_gas (K)
                    1e-16,     // omega0 (rad/s)
                    3.156e16,  // t_age (s, ~1 Gyr)
                    "galaxy", "HST/VLT"}},

    // PAPER_450 — NGC2264 Cone Nebula (wind erosion F_erode=0.05*(t/3Myr))
    {"NGC2264ConeNebula", {"NGC 2264 Cone Nebula (Young Stellar Region)", "HII star-forming region; wind erosion F_erode=0.05*(t/3Myr); Ug3'=G*20Msun/r_star^2; M=100 Msun; Session 120 (PAPER_450)",
                           1.989e32,  // M (kg, ~100 M_sun)
                           3.31e16,   // r (m, ~1.07 pc)
                           1e33,      // L_X (W, HII emission)
                           1e-3,      // B0 (T, HII field)
                           1e-20,     // rho_gas (kg/m^3)
                           1e4,       // T_gas (K)
                           1e-14,     // omega0 (rad/s)
                           3.156e14,  // t_age (s, ~10 Myr)
                           "hii_region", "HST/Spitzer/Chandra"}},

    // PAPER_451 — UGC10214 Tadpole Galaxy (tidal tail v_tail=400 km/s, M_dwarf=3.5e9 Msun)
    {"UGC10214Tadpole", {"UGC 10214 Tadpole Galaxy (tidal encounter)", "Tidal interaction; F_tidal=G*M_dwarf*exp(-t/tau_merge)/d^2; v_tail=400 km/s; tau_merge=2.5e8 yr; Session 120 (PAPER_451)",
                          1.989e41,  // M (kg, ~1e11 M_sun)
                          1.697e21,  // r (m, ~55 kpc)
                          1e37,      // L_X (W)
                          1e-5,      // B0 (T)
                          1e-21,     // rho_gas (kg/m^3)
                          1e4,       // T_gas (K)
                          1e-16,     // omega0 (rad/s)
                          3.156e15,  // t_age (s, ~100 Myr encounter)
                          "galaxy", "HST/ACS"}},

    // PAPER_452 — NGC4676 The Mice (THz H_eff(z)=H(z)*(1+f_THz*log(1+z)), Ug2_THz)
    {"NGC4676TheMice", {"NGC 4676 The Mice (merging galaxy pair)", "Interacting pair; THz Aether expansion H_eff(z)=H(z)*(1+f_THz*log(1+z)); Ug2_THz time-growing; f_THz=0.05; z=0.022; Session 120 (PAPER_452)",
                         1.989e41,  // M (kg, ~1e11 M_sun each; 2x5e10)
                         1.543e21,  // r (m, ~50 kpc)
                         1e37,      // L_X (W)
                         1e-5,      // B0 (T)
                         1e-21,     // rho_gas (kg/m^3)
                         1e4,       // T_gas (K)
                         1e-16,     // omega0 (rad/s)
                         3.156e15,  // t_age (s, ~100 Myr merge stage)
                         "galaxy", "HST/ACS/Spitzer"}},

    // PAPER_453 — NGC6537 Red Spider Nebula (frequency-domain gravity, a=f_total*lambda_P/(2*pi))
    {"NGC6537RedSpider", {"NGC 6537 Red Spider Nebula (bipolar PN)", "Bipolar planetary nebula; freq-domain gravity a=f_total*lambda_Planck/(2*pi); f_super=1.411e16 Hz; t_age=1900 yr; Session 120 (PAPER_453)",
                           1.193e30,  // M (kg, ~0.6 M_sun central WD)
                           7.1e15,    // r (m, bipolar lobe radius)
                           1e36,      // L_X (W, hot PN)
                           1e-4,      // B0 (T)
                           1e-22,     // rho_gas (kg/m^3, bipolar lobe)
                           1e5,       // T_gas (K, fast wind plasma)
                           1e-12,     // omega0 (rad/s)
                           5.996e10,  // t_age (s, ~1900 yr)
                           "planetary_nebula", "HST/NICMOS/VLT"}},

    // PAPER_454 — SMBH Binary LISA Source (Peters inspiral r(t)=r0*(1-t/t_coal)^(1/4))
    {"SMBHBinaryLISA", {"SMBH Binary LISA Gravitational Wave Source", "SMBH binary inspiral; r(t)=r0*(1-t/t_coal)^(1/4) Peters formula; M1=4e6 Msun + M2=2e6 Msun; t_coal=1.555e7 s; freq-domain gravity; Session 120 (PAPER_454)",
                         1.193e37,  // M (kg, ~6e6 M_sun combined)
                         9.461e14,  // r (m, ~0.1 ly initial separation)
                         1e44,      // L_X (W, AGN combined)
                         1e-4,      // B0 (T, accretion disk)
                         1e-20,     // rho_gas (kg/m^3, AGN torus gas)
                         1e8,       // T_gas (K, AGN accretion)
                         1e-4,      // omega0 (rad/s, SMBH orbit)
                         1.555e7,   // t_age (s, coalescence time)
                         "smbh_binary", "LISA/VLBI/HST"}}};

// ============================================================================
// HELPER FUNCTIONS
// ============================================================================

// Get system parameters by name
inline const ObservationalSystem *getSystem(const std::string &name)
{
    auto it = OBSERVATIONAL_SYSTEMS.find(name);
    if (it != OBSERVATIONAL_SYSTEMS.end())
    {
        return &(it->second);
    }
    return nullptr;
}

// Convert system to parameter map for use with PhysicsTerm::compute()
inline std::map<std::string, double> systemToParams(const std::string &system_name)
{
    const ObservationalSystem *sys = getSystem(system_name);
    if (!sys)
        return {};

    std::map<std::string, double> params;
    params["M"] = sys->M;
    params["r"] = sys->r;
    params["L_X"] = sys->L_X;
    params["B0"] = sys->B0;
    params["rho_gas"] = sys->rho_gas;
    params["T_gas"] = sys->T_gas;
    params["omega0"] = sys->omega0;
    params["t_age"] = sys->t_age;

    // Add universal constants
    params["G"] = 6.6743e-11;
    params["c"] = 3e8;
    params["hbar"] = 1.0546e-34;
    params["k_B"] = 1.38e-23;
    params["m_e"] = 9.11e-31;

    return params;
}

// List all available systems
inline std::vector<std::string> listSystems()
{
    std::vector<std::string> names;
    for (const auto &pair : OBSERVATIONAL_SYSTEMS)
    {
        names.push_back(pair.first);
    }
    return names;
}

// Get systems by category
inline std::vector<std::string> getSystemsByCategory(const std::string &category)
{
    std::vector<std::string> names;
    for (const auto &pair : OBSERVATIONAL_SYSTEMS)
    {
        if (pair.second.category == category)
        {
            names.push_back(pair.first);
        }
    }
    return names;
}

#endif // OBSERVATIONAL_SYSTEMS_CONFIG_H
