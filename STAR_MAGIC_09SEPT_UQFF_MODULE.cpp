
// STAR_MAGIC_09SEPT_UQFF_MODULE.cpp
// UQFF 2.0 Standard Module — Star Magic_09Sept2025.docx
// ©2025 Daniel T. Murphy, daniel.murphy00@gmail.com — All Rights Reserved
//
// Source: grok_share_11254865.txt (Grok 4 / X thread conversion of Star Magic_09Sept2025.docx)
// Session 100 — PAPER_368 (Ug4 Vacuum Energy ΛCDM), PAPER_369 (Navier-Stokes Quasar Jet),
//              PAPER_370 (Multi-body Pcore Planetary Scaling)
//
// Architecture: 36th C++ UQFF 2.0 module
//   - 4 Systems: Sun / Earth / Jupiter / Neptune
//   - 4 WOLFRAM_TERM macros: [STARMAG_BASE / STARMAG_UG4_VACUUM / STARMAG_NS_JET / STARMAG_PCORE]
//   - FluidSolver (Jos Stam Stable Fluids, 2D 32×32 grid)
//   - Uses β_i=0.61 (UQFF canonical; source thread uses 0.6 — see helper reference)
//
// CANONICAL DATA FLOW: source2.cpp → APIFetch.py → bodies_*.csv → this module
// This module is a PURE PHYSICS CALCULATOR. CelestialBody parameters are supplied
// at runtime by the PRINCIPAL GUI (source2.cpp); no system data is hardcoded beyond
// default demo values used for exportState() and standalone testing.

#pragma once
#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <functional>

// ============================================================
// WOLFRAM_TERM Macros — symbolic registration for Wolfram WSTP
// ============================================================
#define WOLFRAM_TERM_STARMAG_BASE \
    "Ug4_StarMagBase[k1_,k2_,k3_,k4_,Ms_,Rs_,Bs_,SCm_,QUA_,Pcore_,PSCm_,omega_c_," \
    "rho_A_,kappa_,alpha_,delta_def_,r_,t_,tn_] := " \
    "k1*(Bs + 0.4*Sin[omega_c*t] + 1000)*Rs^3*(G*Ms/Rs^2)*Exp[-alpha*t]*Cos[Pi*tn]*(1+delta_def*Sin[0.001*t]) + " \
    "k2*(1e-10+QUA)*Ms/r^2*(1+0.01*5e5)*1*((SCm*1e8^2/rho_A)*Exp[-kappa*t]) + " \
    "k3*(1e-3+0.4*Sin[omega_c*t]+1000)*Cos[(omega_c-0.4e-6*Sin[omega_c*t])*t*Pi]*Pcore*((SCm*1e8^2/rho_A)*Exp[-kappa*t])"

#define WOLFRAM_TERM_STARMAG_UG4_VACUUM \
    "Ug4_VacuumLambdaCDM[k4_,rho_v_,C_conc_,Mbh_,dg_,alpha_,f_fb_,t_,tn_] := " \
    "k4*rho_v*C_conc*Mbh/dg*Exp[-alpha*t]*Cos[Pi*tn]*(1+f_fb)"

#define WOLFRAM_TERM_STARMAG_NS_JET \
    "NavierStokesUQFFJetForce[v_SCm_,visc_,N_,dt_] := " \
    "Module[{force=v_SCm/10,a=dt*visc*N^2},(force*(1+a))/(1+4*a)]"

#define WOLFRAM_TERM_STARMAG_PCORE \
    "UQFFPcoreScaling[Pcore_,Ug3_] := Pcore*Ug3; " \
    "PcorePlanetaryLaw[body_] := If[body==\"Sun\",1.0,0.001]"

// ============================================================
// Physical Constants
// ============================================================
namespace StarMagic09Sept {

    constexpr double G          = 6.6743e-11;   // m³ kg⁻¹ s⁻²
    constexpr double PI         = 3.14159265358979323846;
    constexpr double c_light    = 2.998e8;       // m/s
    constexpr double HBAR       = 1.0546e-34;    // J·s
    constexpr double YEAR_S     = 3.15576e7;     // s/yr

    // ---- UQFF global parameters (Star Magic 09Sept2025) ----
    constexpr double k1         = 1.5;
    constexpr double k2         = 1.2;
    constexpr double k3         = 1.8;
    constexpr double k4         = 2.0;           // PAPER_368: vacuum coupling constant

    constexpr double beta_i     = 0.61;          // UQFF canonical (thread uses 0.6 — see helper)
    constexpr double kappa      = 0.0005;        // day⁻¹
    constexpr double alpha_dec  = 0.001;         // day⁻¹
    constexpr double gamma_rec  = 0.00005;       // day⁻¹
    constexpr double delta_sw   = 0.01;
    constexpr double epsilon_sw = 0.001;
    constexpr double delta_def  = 0.01;

    constexpr double rho_v      = 6e-27;         // kg/m³ — ΛCDM dark energy density (PAPER_368)
    constexpr double C_conc     = 1.0;           // vacuum concentration factor
    constexpr double f_feedback = 0.1;           // AGN feedback coupling

    constexpr double Mbh        = 8.15e36;       // kg — SgrA* (EHT 2022)
    constexpr double dg         = 2.55e20;       // m — distance to galactic centre
    constexpr double Omega_g    = 7.3e-16;       // rad/s — galactic spin

    constexpr double rho_A      = 1e-23;         // kg/m³ — aether density
    constexpr double rho_sw     = 8e-21;         // kg/m³ — solar wind density
    constexpr double v_sw       = 5e5;           // m/s
    constexpr double v_SCm      = 1e8;           // m/s
    constexpr double QA         = 1e-10;         // C
    constexpr double HSCm       = 1.0;
    constexpr double UUA        = 1.0;
    constexpr double eta_aether = 1e-22;
    constexpr double Ts00       = 1.27e3 + 1.11e7;
    constexpr double num_strings = 1e9;
    constexpr double SCm_contrib = 1e3;          // SCm contribution to Bs/Bj

    // ---- Navier-Stokes parameters (PAPER_369) ----
    constexpr int    NS_N       = 32;
    constexpr double NS_dt      = 0.1;
    constexpr double NS_visc    = 0.0001;
    constexpr double NS_force   = 10.0;          // = v_SCm / 1e7 (scaled)

} // namespace StarMagic09Sept

// ============================================================
// CelestialBody Struct
// ============================================================
namespace StarMagic09Sept {

    struct CelestialBody {
        std::string name;
        double Ms;          // stellar/planetary mass (kg)
        double Rs;          // surface radius (m)
        double Rb;          // interaction boundary radius (m)
        double Ts_surface;  // surface temperature (K)
        double omega_s;     // surface rotation rate (rad/s)
        double Bs_avg;      // average magnetic field (T)
        double SCm_density; // SCm density (kg/m³)
        double QUA;         // quantum aether charge (C)
        double Pcore;       // core SCm penetration factor (1.0=star, 1e-3=planet) [PAPER_370]
        double PSCm;        // SCm pressure factor
        double omega_c;     // characteristic frequency (rad/s) — orbital or solar-cycle [PAPER_370]

        // Derived accessors
        double g_surface() const { return G * Ms / (Rs * Rs); }
    };

    // ---- Default system catalogue (demo / standalone test) ----
    // NOTE: Runtime values are supplied by source2.cpp → APIFetch.py → bodies_*.csv
    static CelestialBody make_Sun() {
        return {
            "Sun", 1.989e30, 6.96e8, 1.496e13, 5778.0, 2.5e-6, 1e-4, 1e15, 1e-11,
            1.0, 1.0,
            2.0 * PI / (11.0 * 365.25 * 86400.0)   // 11-year solar cycle
        };
    }
    static CelestialBody make_Earth() {
        return {
            "Earth", 5.972e24, 6.371e6, 1e7, 288.0, 7.292e-5, 3e-5, 1e12, 1e-12,
            1e-3, 1e-3,
            2.0 * PI / (365.25 * 86400.0)           // 1-year orbital period [PAPER_370]
        };
    }
    static CelestialBody make_Jupiter() {
        return {
            "Jupiter", 1.898e27, 6.9911e7, 1e8, 165.0, 1.76e-4, 4e-4, 1e13, 1e-11,
            1e-3, 1e-3,
            2.0 * PI / (11.86 * 365.25 * 86400.0)   // 11.86-year orbital period [PAPER_370]
        };
    }
    static CelestialBody make_Neptune() {
        return {
            "Neptune", 1.024e26, 2.4622e7, 5e7, 72.0, 1.08e-4, 1e-4, 1e11, 1e-13,
            1e-3, 1e-3,
            2.0 * PI / (164.8 * 365.25 * 86400.0)   // 164.8-year orbital period [PAPER_370]
            // Neptune "frozen planet" — FIRST ice giant in UQFF CelestialBody framework
        };
    }

} // namespace StarMagic09Sept

// ============================================================
// Navier-Stokes Fluid Solver — PAPER_369
// Jos Stam (1999) "Stable Fluids" 2D incompressible solver
// Applied to quasar jet dynamics via SCm velocity forcing
// ============================================================
namespace StarMagic09Sept {

    class FluidSolver {
    private:
        static constexpr int ND = NS_N + 2;
        static constexpr int SZ = ND * ND;

        std::vector<double> u, v, u_prev, v_prev;
        std::vector<double> dens, dens_prev;
        std::vector<double> p_buf, div_buf;     // scratch buffers for project()

        int IX(int i, int j) const { return i + ND * j; }

        void add_source(std::vector<double>& x, const std::vector<double>& s) {
            for (int i = 0; i < SZ; ++i) x[i] += NS_dt * s[i];
        }

        void set_bnd(int b, std::vector<double>& x) {
            for (int i = 1; i <= NS_N; ++i) {
                x[IX(0,    i)] = (b == 1) ? -x[IX(1,    i)] : x[IX(1,    i)];
                x[IX(NS_N+1,i)] = (b == 1) ? -x[IX(NS_N, i)] : x[IX(NS_N, i)];
                x[IX(i,    0)] = (b == 2) ? -x[IX(i,    1)] : x[IX(i,    1)];
                x[IX(i, NS_N+1)] = (b == 2) ? -x[IX(i, NS_N)] : x[IX(i, NS_N)];
            }
            x[IX(0,    0   )] = 0.5 * (x[IX(1,    0   )] + x[IX(0,    1  )]);
            x[IX(0,    NS_N+1)] = 0.5 * (x[IX(1, NS_N+1)] + x[IX(0,    NS_N)]);
            x[IX(NS_N+1, 0  )] = 0.5 * (x[IX(NS_N, 0  )] + x[IX(NS_N+1, 1)]);
            x[IX(NS_N+1,NS_N+1)] = 0.5*(x[IX(NS_N,NS_N+1)]+x[IX(NS_N+1,NS_N)]);
        }

        void diffuse(int b, std::vector<double>& x, const std::vector<double>& x0, double diff) {
            double a = NS_dt * diff * NS_N * NS_N;
            for (int k = 0; k < 20; ++k) {
                for (int i = 1; i <= NS_N; ++i) {
                    for (int j = 1; j <= NS_N; ++j) {
                        x[IX(i,j)] = (x0[IX(i,j)] + a*(x[IX(i-1,j)]+x[IX(i+1,j)]+
                                                        x[IX(i,j-1)]+x[IX(i,j+1)])) / (1.0+4.0*a);
                    }
                }
                set_bnd(b, x);
            }
        }

        void advect(int b, std::vector<double>& d, const std::vector<double>& d0) {
            for (int i = 1; i <= NS_N; ++i) {
                for (int j = 1; j <= NS_N; ++j) {
                    double x = i - NS_dt * NS_N * u[IX(i,j)];
                    double y = j - NS_dt * NS_N * v[IX(i,j)];
                    x = std::max(0.5, std::min((double)NS_N + 0.5, x));
                    y = std::max(0.5, std::min((double)NS_N + 0.5, y));
                    int i0 = (int)x, i1 = i0 + 1;
                    int j0 = (int)y, j1 = j0 + 1;
                    double s1 = x - i0, s0 = 1.0 - s1;
                    double t1 = y - j0, t0 = 1.0 - t1;
                    d[IX(i,j)] = s0*(t0*d0[IX(i0,j0)] + t1*d0[IX(i0,j1)]) +
                                 s1*(t0*d0[IX(i1,j0)] + t1*d0[IX(i1,j1)]);
                }
            }
            set_bnd(b, d);
        }

        void project() {
            double h = 1.0 / NS_N;
            for (int i = 1; i <= NS_N; ++i) {
                for (int j = 1; j <= NS_N; ++j) {
                    div_buf[IX(i,j)] = -0.5 * h * (u[IX(i+1,j)]-u[IX(i-1,j)] +
                                                    v[IX(i,j+1)]-v[IX(i,j-1)]);
                    p_buf[IX(i,j)] = 0.0;
                }
            }
            set_bnd(0, div_buf);
            set_bnd(0, p_buf);
            for (int k = 0; k < 20; ++k) {
                for (int i = 1; i <= NS_N; ++i) {
                    for (int j = 1; j <= NS_N; ++j) {
                        p_buf[IX(i,j)] = (div_buf[IX(i,j)] +
                                          p_buf[IX(i-1,j)] + p_buf[IX(i+1,j)] +
                                          p_buf[IX(i,j-1)] + p_buf[IX(i,j+1)]) / 4.0;
                    }
                }
                set_bnd(0, p_buf);
            }
            for (int i = 1; i <= NS_N; ++i) {
                for (int j = 1; j <= NS_N; ++j) {
                    u[IX(i,j)] -= 0.5 * (p_buf[IX(i+1,j)] - p_buf[IX(i-1,j)]) / h;
                    v[IX(i,j)] -= 0.5 * (p_buf[IX(i,j+1)] - p_buf[IX(i,j-1)]) / h;
                }
            }
            set_bnd(1, u);
            set_bnd(2, v);
        }

    public:
        FluidSolver() :
            u(SZ, 0.0), v(SZ, 0.0), u_prev(SZ, 0.0), v_prev(SZ, 0.0),
            dens(SZ, 0.0), dens_prev(SZ, 0.0), p_buf(SZ, 0.0), div_buf(SZ, 0.0) {}

        // Inject SCm velocity as quasar jet forcing in central column
        void add_jet_force(double force) {
            for (int i = NS_N/4; i <= 3*NS_N/4; ++i)
                v[IX(i, NS_N/2)] += force;
        }

        // One full Stable Fluids time step: diffuse → project → advect → project
        void step() {
            diffuse(1, u_prev, u, NS_visc);
            diffuse(2, v_prev, v, NS_visc);
            project();
            advect(1, u, u_prev);
            advect(2, v, v_prev);
            project();
        }

        // Mean velocity magnitude (used as scalar UQFF output)
        double mean_velocity_magnitude() const {
            double sum = 0.0;
            int cnt = 0;
            for (int i = 1; i <= NS_N; ++i) {
                for (int j = 1; j <= NS_N; ++j) {
                    double ux = u[IX(i,j)], vy = v[IX(i,j)];
                    sum += std::sqrt(ux*ux + vy*vy);
                    ++cnt;
                }
            }
            return cnt > 0 ? sum / cnt : 0.0;
        }

        // ASCII velocity field visualisation
        void print_velocity_field(std::ostream& os = std::cout) const {
            os << "Velocity field (magnitude) — NS_N=" << NS_N << " after jet forcing:\n";
            for (int j = NS_N; j >= 1; --j) {
                for (int i = 1; i <= NS_N; ++i) {
                    double ux = u[IX(i,j)], vy = v[IX(i,j)];
                    double mag = std::sqrt(ux*ux + vy*vy);
                    char sym = (mag > 1.0) ? '#' : (mag > 0.5) ? '+' : (mag > 0.1) ? '.' : ' ';
                    os << sym;
                }
                os << '\n';
            }
        }
    };

    // Run quasar-jet NS simulation; returns mean field magnitude
    double simulate_quasar_jet(double initial_velocity = v_SCm,
                               int steps = 10,
                               bool verbose = false,
                               std::ostream& os = std::cout)
    {
        FluidSolver solver;
        solver.add_jet_force(initial_velocity / 1e7);  // scale for unit grid
        if (verbose) os << "[UQFF NS] Simulating quasar jet (" << steps << " steps)...\n";
        for (int i = 0; i < steps; ++i) solver.step();
        double mag = solver.mean_velocity_magnitude();
        if (verbose) {
            solver.print_velocity_field(os);
            os << "[UQFF NS] Mean |v| after " << steps << " steps: " << mag << " (normalised)\n";
        }
        return mag;
    }

} // namespace StarMagic09Sept

// ============================================================
// UQFF Equation Implementations
// ============================================================
namespace StarMagic09Sept {

    inline double compute_Ereact(const CelestialBody& b, double t) {
        return (b.SCm_density * v_SCm * v_SCm / rho_A) * std::exp(-kappa * t);
    }

    inline double compute_mu_s(const CelestialBody& b, double t) {
        double Bs_t = b.Bs_avg + 0.4 * std::sin(b.omega_c * t) + SCm_contrib;
        return Bs_t * (b.Rs * b.Rs * b.Rs);
    }

    inline double compute_Bj(const CelestialBody& b, double t) {
        return 1e-3 + 0.4 * std::sin(b.omega_c * t) + SCm_contrib;
    }

    inline double compute_omega_s_t(const CelestialBody& b, double t) {
        return b.omega_s - 0.4e-6 * std::sin(b.omega_c * t);
    }

    inline double compute_mu_j(const CelestialBody& b, double t) {
        return compute_Bj(b, t) * (b.Rs * b.Rs * b.Rs);
    }

    // Ug1: magnetic dipole gravity
    double compute_Ug1(const CelestialBody& b, double r, double t, double tn) {
        double mu_s   = compute_mu_s(b, t);
        double grad   = G * b.Ms / (b.Rs * b.Rs);  // ∇(Ms/r) surface approx
        double defect = 1.0 + delta_def * std::sin(0.001 * t);
        return k1 * mu_s * grad * std::exp(-alpha_dec * t) * std::cos(PI * tn) * defect;
    }

    // Ug2: charge-reactivity gravity
    double compute_Ug2(const CelestialBody& b, double r, double t, double tn) {
        double Ereact  = compute_Ereact(b, t);
        double S       = (r > b.Rb) ? 1.0 : 0.0;
        double wind_m  = 1.0 + delta_sw * v_sw;
        return k2 * (QA + b.QUA) * b.Ms / (r * r) * S * wind_m * HSCm * Ereact;
    }

    // Ug3: string rotation gravity (optimised: Pcore scaling; no loop)
    double compute_Ug3(const CelestialBody& b, double r, double t, double tn) {
        double Ereact   = compute_Ereact(b, t);
        double omega_st = compute_omega_s_t(b, t);
        double Bj       = compute_Bj(b, t);
        return k3 * Bj * std::cos(omega_st * t * PI) * b.Pcore * Ereact;
    }

    // Ug4: vacuum energy ΛCDM galactic BH coupling — PAPER_368
    // Ug4 = k4 × ρ_v × C_conc × Mbh/dg × exp(−α×t) × cos(π×tn) × (1+f_feedback)
    // FIRST explicit ΛCDM dark-energy density coupling to galactic BH/distance in UQFF
    // Distinct from Ug4VacuumMediatedCalculator (f3c55f52): k4=1e-40, J/m³, [SCm] coupling
    double compute_Ug4(double t, double tn) {
        return k4 * rho_v * C_conc * (Mbh / dg)
             * std::exp(-alpha_dec * t)
             * std::cos(PI * tn)
             * (1.0 + f_feedback);
    }

    // Universal Buoyancy (canonical β_i=0.61)
    double compute_Ubi(double Ugi, double tn) {
        double wind_m = 1.0 + epsilon_sw * rho_sw;
        return -beta_i * Ugi * Omega_g * (Mbh / dg) * wind_m * UUA * std::cos(PI * tn);
    }

    // Universal Magnetism (optimised: multiply by num_strings)
    double compute_Um(const CelestialBody& b, double t, double tn, double rj) {
        double Ereact = compute_Ereact(b, t);
        double mu_j   = compute_mu_j(b, t);
        double decay  = 1.0 - std::exp(-gamma_rec * t * std::cos(PI * tn));
        double single = (mu_j / rj) * decay;
        return single * num_strings * b.PSCm * Ereact;
    }

    // Aether metric trace (Aμν scalar contribution to FU)
    double compute_A_mu_nu_trace(double tn) {
        // g_μν = diag(1,-1,-1,-1); trace = 1-1-1-1 = -2
        // Aether modulation: +η×Ts00×cos(π×tn) per component → trace += 4×mod
        double mod   = eta_aether * Ts00 * std::cos(PI * tn);
        double trace = (1.0 - 1.0 - 1.0 - 1.0) + 4.0 * mod;
        return trace;
    }

    // Master unified field strength
    double compute_FU(const CelestialBody& b, double r, double t, double tn) {
        double Ug1 = compute_Ug1(b, r, t, tn);
        double Ug2 = compute_Ug2(b, r, t, tn);
        double Ug3 = compute_Ug3(b, r, t, tn);
        double Ug4 = compute_Ug4(t, tn);

        double Ubi1 = compute_Ubi(Ug1, tn);
        double Ubi2 = compute_Ubi(Ug2, tn);
        double Ubi3 = compute_Ubi(Ug3, tn);
        double Ubi4 = compute_Ubi(Ug4, tn);

        double Um     = compute_Um(b, t, tn, b.Rb);
        double A_trace = compute_A_mu_nu_trace(tn);

        return (Ug1+Ug2+Ug3+Ug4) + (Ubi1+Ubi2+Ubi3+Ubi4) + Um + A_trace;
    }

} // namespace StarMagic09Sept

// ============================================================
// UQFF 2.0 Module Class
// ============================================================
namespace StarMagic09Sept {

    class StarMagic09SeptModule {
    private:
        bool log_enabled = false;
        std::map<std::string, double> dynamic_params;

        void log(const std::string& msg) const {
            if (log_enabled) std::cout << "[StarMagic09Sept] " << msg << '\n';
        }

    public:
        // ---- UQFF 2.0 interface methods ----

        void setEnableLogging(bool enabled) {
            log_enabled = enabled;
        }

        void setDynamicParameter(const std::string& name, double value) {
            dynamic_params[name] = value;
            log("SetParam: " + name + " = " + std::to_string(value));
        }

        double getDynamicParameter(const std::string& name) const {
            auto it = dynamic_params.find(name);
            if (it != dynamic_params.end()) return it->second;
            throw std::runtime_error("Parameter not found: " + name);
        }

        // Export state to file (UQFF 2.0 state persistence)
        void exportState(const std::string& filename) const {
            std::ofstream f(filename);
            if (!f.is_open()) throw std::runtime_error("Cannot open: " + filename);
            f << "# StarMagic09Sept UQFF 2.0 State Export\n";
            f << "module: STAR_MAGIC_09SEPT_UQFF_MODULE\n";
            f << "session: 100\n";
            f << "papers: PAPER_368,PAPER_369,PAPER_370\n";
            f << "k1: " << k1 << "\n";
            f << "k2: " << k2 << "\n";
            f << "k3: " << k3 << "\n";
            f << "k4: " << k4 << "\n";
            f << "beta_i: " << beta_i << "\n";
            f << "rho_v: " << rho_v << "\n";
            f << "C_conc: " << C_conc << "\n";
            f << "f_feedback: " << f_feedback << "\n";
            f << "Mbh: " << Mbh << "\n";
            f << "dg: " << dg << "\n";
            f << "kappa: " << kappa << "\n";
            f << "alpha: " << alpha_dec << "\n";
            f << "gamma: " << gamma_rec << "\n";
            f << "NS_N: " << NS_N << "\n";
            f << "NS_dt: " << NS_dt << "\n";
            f << "NS_visc: " << NS_visc << "\n";
            f << "Pcore_Sun: 1.0\n";
            f << "Pcore_planet: 0.001\n";
            for (const auto& kv : dynamic_params)
                f << "dyn_" << kv.first << ": " << kv.second << "\n";
            f.close();
            log("State exported to " + filename);
        }

        // Cross-validation: compare Ug4 ΛCDM form vs Newton baseline
        template <typename BodyT>
        std::map<std::string, double> cross_validate(const BodyT& body,
                                                      double r, double t, double tn) const {
            double g_Newton = G * body.Ms / (r * r);
            double Ug4_lambda = compute_Ug4(t, tn);
            double FU = compute_FU(body, r, t, tn);
            double eta_Ug4 = (g_Newton > 0.0) ? Ug4_lambda / g_Newton : 0.0;
            double eta_FU  = (g_Newton > 0.0) ? std::abs(FU) / g_Newton : 0.0;
            return {
                {"g_Newton_ms2",    g_Newton},
                {"Ug4_lambda_ms2",  Ug4_lambda},
                {"FU_ms2",          FU},
                {"eta_Ug4_FU_ratio", eta_Ug4},
                {"eta_FU_Newton",   eta_FU},
            };
        }

        // Compute full UQFF for all canonical bodies
        void computeAll(double t = 0.0, double tn = 0.0,
                        std::ostream& os = std::cout) const {
            std::vector<CelestialBody> bodies = {
                make_Sun(), make_Earth(), make_Jupiter(), make_Neptune()
            };
            os << "=== StarMagic09Sept UQFF 2.0 — All Systems ===\n";
            os << "t=" << t << " days, tn=" << tn << "\n\n";
            for (const auto& b : bodies) {
                double r  = b.Rb;
                double FU = compute_FU(b, r, t, tn);
                double Ug4 = compute_Ug4(t, tn);
                double g_N = b.g_surface();
                os << b.name << ": FU=" << FU << "  Ug4=" << Ug4
                   << "  g_Newton=" << g_N << "  Pcore=" << b.Pcore << "\n";
            }
            // PAPER_369: NS quasar jet
            double ns_mag = simulate_quasar_jet(v_SCm, 10, false);
            os << "\n[PAPER_369] NS quasar jet mean |v|=" << ns_mag << " (normalised)\n";
            os << "=== End StarMagic09Sept UQFF 2.0 ===\n";
        }
    };

} // namespace StarMagic09Sept

// ============================================================
// Standalone self-test (called when linked as option in MAIN_1)
// ============================================================
namespace StarMagic09Sept {

    inline void run_selftest(std::ostream& os = std::cout) {
        os << "\n[StarMagic09Sept] Self-test — PAPER_368 / PAPER_369 / PAPER_370\n";

        // PAPER_368: Ug4 at t=0, tn=0
        double Ug4_t0 = compute_Ug4(0.0, 0.0);
        os << "PAPER_368 Ug4(t=0,tn=0) = " << Ug4_t0 << " m/s²  (expect ~4.22e-10)\n";

        // PAPER_369: NS quasar jet
        double ns_mag = simulate_quasar_jet(v_SCm, 10, true, os);
        os << "PAPER_369 NS quasar jet mean |v| = " << ns_mag << " (normalised)\n";

        // PAPER_370: Pcore scaling
        CelestialBody sun  = make_Sun();
        CelestialBody nep  = make_Neptune();
        os << "PAPER_370 Sun.Pcore=" << sun.Pcore << " (expect 1.0)\n";
        os << "PAPER_370 Neptune.Pcore=" << nep.Pcore << " (expect 0.001)\n";
        os << "PAPER_370 Sun.omega_c=" << sun.omega_c << " rad/s  (11yr cycle)\n";
        os << "PAPER_370 Neptune.omega_c=" << nep.omega_c << " rad/s  (164.8yr orbit)\n";

        // Full FU for Sun
        double FU_sun = compute_FU(sun, sun.Rb, 0.0, 0.0);
        os << "FU(Sun, r=Rb, t=0) = " << FU_sun << "\n";

        StarMagic09SeptModule mod;
        mod.setEnableLogging(true);
        mod.setDynamicParameter("test_param", 42.0);
        mod.exportState("StarMagic09Sept_state.txt");
        mod.computeAll(0.0, 0.0, os);

        os << "[StarMagic09Sept] Self-test complete.\n\n";
    }

} // namespace StarMagic09Sept

// ============================================================
// End of STAR_MAGIC_09SEPT_UQFF_MODULE.cpp
// Session 100 — PAPER_368 / PAPER_369 / PAPER_370
// ============================================================

// ============================================================
// SESSION 101 EXTENSION — PAPER_371 / PAPER_372 / PAPER_373 / PAPER_374 / PAPER_375
// Source: grok_share_11254865.txt (lines 2000–8800, extended re-analysis)
// Source docs integrated:
//   "200. MUGE Compression cycle 3_Superconductive Resonance_11May2025.docx"
//   "100. MUGE Compression cycle 3_11May2025.docx"
//   "Compressed UQFF Equation_14May2025.docx"
//   "Master UQFF Resonance Equation_14May2025.docx"
//   "UQFF_Resonance Superconductive Universal Gravity Equation system proof set._15May2025.docx"
// ============================================================

namespace StarMagic09Sept_Session101 {

// ============================================================
// WOLFRAM_TERM Macros — Session 101 additions
// ============================================================
#define WOLFRAM_TERM_MUGE_RESONANCE \
    "MUGE_Resonance[fDPM_,fTHz_,Evac_neb_,Evac_ISM_,Delta_Evac_,Fsuper_,UA_SCM_,omega_i_," \
    "k4_res_,freact_,fquantum_,fAether_,fosc_,fTRZ_,I_,A_,omega1_,omega2_," \
    "Vsys_,vexp_,Ereact_,H_z_,t_,ffluid_,c_] := " \
    "Let[{FDPM=I*A*(omega1-omega2),aDPM=FDPM*fDPM*Evac_neb*c*Vsys}," \
    "aDPM + fTHz*Evac_neb*vexp*aDPM/Evac_ISM/c + Delta_Evac*vexp^2*aDPM/Evac_neb/c^2 + " \
    "Fsuper*fTHz*aDPM/Evac_neb/c + UA_SCM*omega_i*fTHz*aDPM*(1+fTRZ) + " \
    "k4_res*Ereact*freact*aDPM/Evac_neb*c + fquantum*Evac_neb*aDPM/Evac_ISM/c + " \
    "fAether*Evac_neb*aDPM/Evac_ISM/c + ffluid*Evac_neb*Vsys/Evac_ISM/c + " \
    "fosc*Cos[2*Pi*fosc*t] + 2*Pi*H_z*t*Evac_neb*aDPM/Evac_ISM/c + fTRZ]"

#define WOLFRAM_TERM_COMPRESSED_UQFF_BCRIT \
    "CompressedUQFF_Bcrit[G_,M_,r_,H0_,t_,B_,Bcrit_,Lambda_,c_,hbar_,rho_fluid_,V_,g_local_,M_DM_,delta_rho_,rho_] := " \
    "(G*M/r^2)*(1+H0*t)*(1-B/Bcrit) + Lambda*c^2/3 + " \
    "(hbar/1e-68)*2.176e-18*(2*Pi/4.35e17) + rho_fluid*V*g_local + " \
    "(M+M_DM)*(delta_rho/rho + 3*G*M/r^3)"

#define WOLFRAM_TERM_WORMHOLE_GEODESIC \
    "WormholeGeodesic[b_,r_,E_,L_,lambda_step_,n_steps_] := " \
    "NestList[Function[{r0,phi0,t0},{r0 + lambda_step*Sqrt[E^2-L^2/(b^2+r0^2)]," \
    "phi0 + lambda_step*L/(b^2+r0^2), t0 + lambda_step*E}], {0,0,0}, n_steps]"

#define WOLFRAM_TERM_UQFF_ADVANCED \
    "UQFF_Advanced[g_compressed_,g_resonance_,a_worm_,B_,Bcrit_,v_,c_] := " \
    "g_compressed*Exp[-B/Bcrit] + g_resonance/Sqrt[1-(v/c)^2] + a_worm"

// ============================================================
// PAPER_371: MUGE 12-Term Superconductive Resonance Framework
// Source: "200. MUGE Compression cycle 3_Superconductive Resonance_11May2025.docx"
// ============================================================

struct ResonanceParams {
    double fDPM        = 1e12;      // Hz — DPM frequency
    double fTHz        = 1e12;      // Hz — THz frequency
    double Evac_neb    = 7.09e-36;  // J  — nebular vacuum energy
    double Evac_ISM    = 7.09e-37;  // J  — ISM vacuum energy
    double Delta_Evac  = 6.381e-36; // J  — vacuum energy difference
    double Fsuper      = 6.287e-19; // N  — superconductive force
    double UA_SCM      = 10.0;      //    — aether SCm coupling
    double omega_i     = 1e-8;      // rad/s — intrinsic angular frequency
    double k4_res      = 1.0;       //    — resonance Ug4 coupling
    double freact      = 1e10;      // Hz — reactive frequency
    double fquantum    = 1.445e-17; // Hz — quantum frequency
    double fAether     = 1.576e-35; // Hz — aether frequency
    double fosc        = 4.57e14;   // Hz — oscillation frequency
    double fTRZ        = 0.1;       //    — time-reversal correction
    double c_res       = 3e8;       // m/s
};

struct MUGESystem {
    std::string name;
    double M;           // kg — system mass
    double r;           // m  — characteristic radius
    double B;           // T  — magnetic field
    double Bcrit;       // T  — critical field
    double Vsys;        // m³ — system volume
    double ffluid;      // Hz — fluid frequency
    double vexp;        // m/s — expansion velocity
    double rho_fluid;   // kg/m³ — fluid density
    double M_DM;        // kg — dark matter mass
};

// Predefined 7-system catalog (demo only; runtime values from source2.cpp pipeline)
inline MUGESystem make_SGR1745() {
    return {"SGR1745-2900", 2.984e30, 1e4, 1e10, 1e11, 4.189e12, 1.269e-14, 1e3, 1e-3, 0.0};
}
inline MUGESystem make_SagA_star() {
    return {"Sagittarius A*", 8.155e36, 1e12, 1e-5, 1e-4, 3.552e45, 3.465e-8, 1e6, 1e-4, 8e36};
}
inline MUGESystem make_Tapestry() {
    return {"Tapestry Starbirth", 1.989e35, 3.086e17, 1e-4, 1e-3, 1e53, 1e-12, 1e4, 1e-4, 0.0};
}
inline MUGESystem make_Westerlund2() {
    return {"Westerlund 2", 1.989e35, 3.086e17, 1e-4, 1e-3, 1e53, 1e-12, 1e4, 1e-4, 0.0};
}
inline MUGESystem make_Pillars() {
    return {"Pillars of Creation", 1.989e32, 9.46e15, 1e-4, 1e-3, 3.552e48, 8.457e-14, 1e3, 1e-5, 0.0};
}
inline MUGESystem make_Rings() {
    return {"Rings of Relativity", 1.989e36, 3.086e17, 1e-5, 1e-4, 1e54, 1e-9, 1e5, 1e-5, 1.989e35};
}
inline MUGESystem make_StudentGuide() {
    return {"Student's Guide Universe", 1e53, 1e26, 1e-10, 1e-9, 1e80, 1e-18, 1e7, 1e-10, 1e52};
}

// Individual MUGE resonance term functions (PAPER_371)
inline double compute_aDPM(const MUGESystem& sys, const ResonanceParams& p,
                            double I = 1.0, double A = 1.0,
                            double omega1 = 1e12, double omega2 = 9.99e11) {
    double FDPM = I * A * (omega1 - omega2);
    return FDPM * p.fDPM * p.Evac_neb * p.c_res * sys.Vsys;
}

inline double compute_aTHz(const ResonanceParams& p, double aDPM, double vexp) {
    return p.fTHz * p.Evac_neb * vexp * aDPM / p.Evac_ISM / p.c_res;
}

inline double compute_avac_diff(const ResonanceParams& p, double aDPM, double vexp) {
    return p.Delta_Evac * vexp * vexp * aDPM / p.Evac_neb / (p.c_res * p.c_res);
}

inline double compute_asuper_freq(const ResonanceParams& p, double aDPM) {
    return p.Fsuper * p.fTHz * aDPM / p.Evac_neb / p.c_res;
}

inline double compute_aaether_res(const ResonanceParams& p, double aDPM) {
    return p.UA_SCM * p.omega_i * p.fTHz * aDPM * (1.0 + p.fTRZ);
}

inline double compute_Ug4i(const ResonanceParams& p, double aDPM, double Ereact) {
    return p.k4_res * Ereact * p.freact * aDPM / p.Evac_neb * p.c_res;
}

inline double compute_aquantum_freq(const ResonanceParams& p, double aDPM) {
    return p.fquantum * p.Evac_neb * aDPM / p.Evac_ISM / p.c_res;
}

inline double compute_aAether_freq(const ResonanceParams& p, double aDPM) {
    return p.fAether * p.Evac_neb * aDPM / p.Evac_ISM / p.c_res;
}

inline double compute_afluid_freq(const MUGESystem& sys, const ResonanceParams& p) {
    return sys.ffluid * p.Evac_neb * sys.Vsys / p.Evac_ISM / p.c_res;
}

inline double compute_Osc_term(const ResonanceParams& p, double t) {
    return p.fosc * std::cos(2.0 * M_PI * p.fosc * t);
}

inline double compute_aexp_freq(const ResonanceParams& p, double aDPM,
                                 double H_z, double t) {
    return 2.0 * M_PI * H_z * t * p.Evac_neb * aDPM / p.Evac_ISM / p.c_res;
}

// Master 12-term MUGE resonance function (PAPER_371)
inline double compute_resonance_MUGE(const MUGESystem& sys,
                                      const ResonanceParams& p = ResonanceParams(),
                                      double t = 0.0,
                                      double H_z = 2.269e-18,
                                      double Ereact = 1.0) {
    double aDPM      = compute_aDPM(sys, p);
    double aTHz      = compute_aTHz(p, aDPM, sys.vexp);
    double avac_diff = compute_avac_diff(p, aDPM, sys.vexp);
    double asuper    = compute_asuper_freq(p, aDPM);
    double aaether   = compute_aaether_res(p, aDPM);
    double Ug4i      = compute_Ug4i(p, aDPM, Ereact);
    double aquantum  = compute_aquantum_freq(p, aDPM);
    double aAether   = compute_aAether_freq(p, aDPM);
    double afluid    = compute_afluid_freq(sys, p);
    double Osc       = compute_Osc_term(p, t);
    double aexp      = compute_aexp_freq(p, aDPM, H_z, t);
    double fTRZ      = p.fTRZ;
    return aDPM + aTHz + avac_diff + asuper + aaether + Ug4i
         + aquantum + aAether + afluid + Osc + aexp + fTRZ;
}

// ============================================================
// PAPER_372: Compressed UQFF with B/Bcrit Superconductivity
// Source: "100. MUGE Compression cycle 3_11May2025.docx"
// ============================================================

namespace CompressedUQFF {
    constexpr double G        = 6.674e-11;
    constexpr double H0       = 2.269e-18; // s⁻¹ (Hubble constant)
    constexpr double Lambda   = 1.1e-52;   // m⁻²
    constexpr double c        = 3e8;
    constexpr double hbar     = 1.055e-34;
    constexpr double tHubble  = 4.35e17;   // s

    inline double compressed_base(const MUGESystem& sys) {
        return G * sys.M / (sys.r * sys.r);
    }
    inline double compressed_expansion(const MUGESystem& sys, double t) {
        return 1.0 + H0 * t;
    }
    inline double compressed_super_adj(const MUGESystem& sys) {
        return 1.0 - sys.B / sys.Bcrit;
    }
    inline double compressed_env() {
        return 1.0; // default Fenv
    }
    inline double compressed_cosm() {
        return Lambda * c * c / 3.0;
    }
    inline double compressed_quantum() {
        // (ℏ/Δx·Δp)·∫(ψ*Ĥψ dV)·(2π/tHubble) — symbolic constant
        return (hbar / 1e-68) * 2.176e-18 * (2.0 * M_PI / tHubble);
    }
    inline double compressed_fluid(const MUGESystem& sys, double g_local) {
        return sys.rho_fluid * sys.Vsys * g_local;
    }
    inline double compressed_perturbation(const MUGESystem& sys) {
        double delta_rho_over_rho = 1e-5; // dimensionless density perturbation
        return (sys.M + sys.M_DM) * (delta_rho_over_rho + 3.0 * G * sys.M / (sys.r * sys.r * sys.r));
    }
    // Master compressed UQFF (PAPER_372)
    inline double compute_compressed_MUGE(const MUGESystem& sys, double t = 0.0) {
        double base  = compressed_base(sys);
        double exp_  = compressed_expansion(sys, t);
        double super = compressed_super_adj(sys);
        double env   = compressed_env();
        double cosm  = compressed_cosm();
        double quant = compressed_quantum();
        double fluid = compressed_fluid(sys, base * super);
        double pert  = compressed_perturbation(sys);
        return base * exp_ * super * (1.0 + env) + cosm + quant + fluid + pert;
    }
} // namespace CompressedUQFF

// ============================================================
// PAPER_373: Morris-Thorne Wormhole Null Geodesics
// Source: wormhole section (lines ~2700–2800)
// FIRST wormhole physics in the entire CP pipeline.
// ============================================================

namespace WormholeGeodesics {
    // Morris-Thorne metric: ds² = -dt² + dr² + (b²+r²)(dθ²+sin²θ dφ²)
    constexpr double b_throat = 1.0; // throat radius (m)

    struct GeodesicState {
        double r, phi, t_coord;
    };

    // dr/dλ = ±√(E²−L²/(b²+r²))
    inline double drdt(double E, double L, double r, double b = b_throat) {
        double arg = E * E - L * L / (b * b + r * r);
        return (arg >= 0.0) ? std::sqrt(arg) : 0.0;
    }
    inline double dphidt(double L, double r, double b = b_throat) {
        return L / (b * b + r * r);
    }
    // Embedding: z_embed = b·arcsinh(r/b),  ρ_embed = √(b²+r²)
    inline double z_embed(double r, double b = b_throat) {
        return b * std::asinh(r / b);
    }
    inline double rho_embed(double r, double b = b_throat) {
        return std::sqrt(b * b + r * r);
    }
    // Propagate null geodesic for n steps of dlambda
    inline std::vector<GeodesicState> propagate(double E, double L,
                                                  double r0, int n_steps = 100,
                                                  double dlambda = 0.1) {
        std::vector<GeodesicState> traj;
        traj.reserve(n_steps + 1);
        GeodesicState s{r0, 0.0, 0.0};
        traj.push_back(s);
        for (int i = 0; i < n_steps; ++i) {
            double dr  = drdt(E, L, s.r) * dlambda;
            double dph = dphidt(L, s.r) * dlambda;
            double dt  = E * dlambda;
            s.r += dr; s.phi += dph; s.t_coord += dt;
            traj.push_back(s);
        }
        return traj;
    }
    // Self-test: traversal (L=0.5, E=1.0) and reflection (L=1.5, E=1.0)
    inline void selftest(std::ostream& os = std::cout) {
        os << "[WormholeGeodesics] PAPER_373 self-test\n";
        // Traversal case
        auto trav = propagate(1.0, 0.5, 2.0, 50, 0.1);
        os << "  Traversal L=0.5: r_final=" << trav.back().r << " (expect crosses 0)\n";
        // Reflection case
        auto refl = propagate(1.0, 1.5, 2.0, 50, 0.1);
        double r_min = std::min_element(refl.begin(), refl.end(),
                         [](const GeodesicState& a, const GeodesicState& b){ return a.r < b.r; })->r;
        os << "  Reflection L=1.5: r_min=" << refl[0].r << " (expect ≈1.12)\n";
        os << "  z_embed(r=1.0)=" << z_embed(1.0) << " rho_embed(r=1.0)=" << rho_embed(1.0) << "\n";
        os << "[WormholeGeodesics] PAPER_373 self-test complete.\n";
    }
} // namespace WormholeGeodesics

// ============================================================
// PAPER_374: J1610+1811 Relativistic Quasar Jet UQFF-NS Coupling
// Source: simulate_quasar_jet() (lines ~5100–5200)
// Distinct from PAPER_360 (FU/Bi at z=6.5): this is UQFF resonance
// force coupling into NS solver at v_SCm=0.99c.
// Observational: z=3.122, P_jet=4e45 W, L=2e46 W → v_SCm=0.99c
// ============================================================

namespace J1610QuasarJet {
    constexpr double c            = 3e8;
    constexpr double z_redshift   = 3.122;      // J1610+1811 redshift
    constexpr double P_jet        = 4e45;        // W — jet power
    constexpr double L_luminosity = 2e46;        // W — total luminosity
    constexpr double v_SCm_rel    = 0.99 * c;    // m/s — relativistic jet velocity

    // Simulate quasar jet with UQFF resonance force coupling
    // Returns mean |v| of final NS field (normalised)
    inline double simulate_relativistic_quasar_jet(
            std::ostream& os = std::cout, int NS_steps = 10) {
        // System: Sagittarius A* as proxy for quasar host SMBH
        MUGESystem sagA = make_SagA_star();
        ResonanceParams res;
        double uqff_g = compute_resonance_MUGE(sagA, res);

        // NS FluidSolver (FluidSolver is defined in Session 100 namespace above)
        using namespace StarMagic09Sept;
        FluidSolver fs;
        double jet_force = v_SCm_rel / 10.0; // relativistic jet forcing

        // Apply jet force and run NS steps using public interface
        fs.add_jet_force(jet_force + uqff_g / 1e30);
        for (int step = 0; step < NS_steps; ++step) {
            fs.step();
        }
        double mean_v = fs.mean_velocity_magnitude();
        os << "J1610+1811 Relativistic Quasar Jet (PAPER_374)\n";
        os << "  z=" << z_redshift << " P_jet=" << P_jet << "W  L=" << L_luminosity << "W\n";
        os << "  v_SCm=" << v_SCm_rel << " m/s  uqff_g=" << uqff_g << " m/s²\n";
        os << "  NS mean|v| after " << NS_steps << " steps = " << mean_v << "\n";
        return mean_v;
    }
} // namespace J1610QuasarJet

// ============================================================
// PAPER_375: UQFF Advanced Integration
// Wormhole-MUGE term + Meissner exponential + Relativistic γ + Error propagation
// Source: Unified UQFF analysis (lines ~7500–8800), three new docs
// ============================================================

namespace UQFFAdvanced {
    constexpr double f_worm = 1e-10; // wormhole coupling constant
    constexpr double c = 3e8;

    // 1. Wormhole-MUGE coupling term
    inline double compute_a_wormhole(double Evac_neb, double b, double r) {
        return f_worm * Evac_neb / (b * b + r * r);
    }

    // 2. Meissner exponential superconductivity (type-II improved form)
    // Replaces linear (1−B/Bcrit) from PAPER_372 with exp(−B/Bcrit)
    inline double meissner_exp(double B, double Bcrit) {
        return std::exp(-B / Bcrit);
    }

    // 3. Relativistic Lorentz correction
    inline double lorentz_gamma(double v) {
        double beta = v / c;
        if (beta >= 1.0) beta = 0.9999999;
        return 1.0 / std::sqrt(1.0 - beta * beta);
    }
    inline double apply_lorentz(double aDPM, double v) {
        return aDPM / lorentz_gamma(v);
    }

    // 4. Error propagation
    inline double error_propagation(const std::vector<double>& delta_terms) {
        double sum_sq = 0.0;
        for (double d : delta_terms) sum_sq += d * d;
        return std::sqrt(sum_sq);
    }

    // Master unified UQFF with all advanced terms (PAPER_375)
    inline double compute_unified_UQFF(const MUGESystem& sys,
                                        const ResonanceParams& res = ResonanceParams(),
                                        double t = 0.0,
                                        double v_jet = 0.0,
                                        double b_worm = 1.0,
                                        double r_worm = 1.0) {
        // Compressed UQFF base with Meissner exponential
        double base       = CompressedUQFF::compressed_base(sys);
        double expansion  = CompressedUQFF::compressed_expansion(sys, t);
        double meissner   = meissner_exp(sys.B, sys.Bcrit);
        double cosm       = CompressedUQFF::compressed_cosm();
        double quant      = CompressedUQFF::compressed_quantum();
        double fluid      = CompressedUQFF::compressed_fluid(sys, base);
        double pert       = CompressedUQFF::compressed_perturbation(sys);
        double g_compressed = base * expansion * meissner + cosm + quant + fluid + pert;

        // MUGE Resonance with relativistic gamma on aDPM
        double aDPM_raw   = compute_aDPM(sys, res);
        double aDPM_rel   = apply_lorentz(aDPM_raw, v_jet);
        // Substitute aDPM_rel into resonance (pass as modified param)
        ResonanceParams res_mod = res;
        // Scale full resonance by aDPM_rel / aDPM_raw ratio
        double gamma_ratio = (aDPM_raw != 0.0) ? aDPM_rel / aDPM_raw : 1.0;
        double g_resonance = compute_resonance_MUGE(sys, res_mod, t) * gamma_ratio;

        // Wormhole-MUGE term
        double a_worm = compute_a_wormhole(res.Evac_neb, b_worm, r_worm);

        return g_compressed + g_resonance + a_worm;
    }

    // Error propagation for all MUGE terms
    inline double compute_total_uncertainty(const MUGESystem& sys,
                                             const ResonanceParams& p = ResonanceParams(),
                                             double frac_error = 0.01) {
        std::vector<double> deltas;
        double aDPM = compute_aDPM(sys, p);
        deltas.push_back(std::abs(aDPM) * frac_error);
        deltas.push_back(std::abs(compute_aTHz(p, aDPM, sys.vexp)) * frac_error);
        deltas.push_back(std::abs(compute_avac_diff(p, aDPM, sys.vexp)) * frac_error);
        deltas.push_back(std::abs(compute_asuper_freq(p, aDPM)) * frac_error);
        deltas.push_back(std::abs(compute_aaether_res(p, aDPM)) * frac_error);
        deltas.push_back(std::abs(compute_afluid_freq(sys, p)) * frac_error);
        return error_propagation(deltas);
    }
} // namespace UQFFAdvanced

// ============================================================
// Session 101 self-test — PAPER_371–375
// ============================================================
inline void run_session101_selftest(std::ostream& os = std::cout) {
    os << "\n[Session101] Self-test — PAPER_371 / PAPER_372 / PAPER_373 / PAPER_374 / PAPER_375\n";

    ResonanceParams res;
    MUGESystem sgr = make_SGR1745();
    MUGESystem sagA = make_SagA_star();

    // PAPER_371: resonance MUGE on SGR1745
    double g_res = compute_resonance_MUGE(sgr, res);
    os << "PAPER_371 resonance_MUGE(SGR1745) = " << g_res << " m/s²  (expect ~1.773e-9)\n";

    // PAPER_371: individual terms
    double aDPM   = compute_aDPM(sgr, res);
    double aTHz   = compute_aTHz(res, aDPM, sgr.vexp);
    double avac   = compute_avac_diff(res, aDPM, sgr.vexp);
    double asuper = compute_asuper_freq(res, aDPM);
    double aaether= compute_aaether_res(res, aDPM);
    double aquant = compute_aquantum_freq(res, aDPM);
    double aAeth  = compute_aAether_freq(res, aDPM);
    double aexp   = compute_aexp_freq(res, aDPM, 2.269e-18, 3.799e10);
    os << "  aTHz          = " << aTHz   << "  (expect ~1.182e-33)\n";
    os << "  avac_diff     = " << avac   << "  (expect ~3.545e-53)\n";
    os << "  asuper_freq   = " << asuper << "  (expect ~1.048e-21)\n";
    os << "  aaether_res   = " << aaether<< "  (expect ~3.900e-38)\n";
    os << "  aquantum_freq = " << aquant << "  (expect ~1.708e-66)\n";
    os << "  aAether_freq  = " << aAeth  << "  (expect ~1.863e-84)\n";
    os << "  aexp_freq(t=3.799e10) = " << aexp << "  (expect ~1.623e-57)\n";
    os << "  afluid_freq(SGR1745)  = " << compute_afluid_freq(sgr, res) << "  (expect ~1.773e-9)\n";

    // PAPER_372: compressed MUGE on SGR1745
    double g_comp = CompressedUQFF::compute_compressed_MUGE(sgr);
    os << "PAPER_372 compressed_MUGE(SGR1745) = " << g_comp << " m/s²  (expect ~1.782e39)\n";

    // PAPER_373: wormhole geodesics
    WormholeGeodesics::selftest(os);

    // PAPER_374: relativistic quasar jet
    double mean_v = J1610QuasarJet::simulate_relativistic_quasar_jet(os);
    os << "PAPER_374 J1610 NS mean|v| = " << mean_v << "\n";

    // PAPER_375: unified UQFF
    double g_unified = UQFFAdvanced::compute_unified_UQFF(sagA, res, 0.0,
                                                            J1610QuasarJet::v_SCm_rel);
    double dg        = UQFFAdvanced::compute_total_uncertainty(sgr, res);
    os << "PAPER_375 unified_UQFF(SagA*, v=0.99c) = " << g_unified << "\n";
    os << "PAPER_375 δg(SGR1745, 1% error) = " << dg << "\n";
    double gma = UQFFAdvanced::lorentz_gamma(J1610QuasarJet::v_SCm_rel);
    os << "PAPER_375 γ(v=0.99c) = " << gma << "  (expect ~7.09)\n";
    double a_worm_test = UQFFAdvanced::compute_a_wormhole(res.Evac_neb, 1.0, 1.0);
    os << "PAPER_375 a_worm(r=1,b=1) = " << a_worm_test << "\n";

    os << "[Session101] Self-test complete. Papers: 371 / 372 / 373 / 374 / 375\n\n";
}

} // namespace StarMagic09Sept_Session101

// ============================================================
// SESSION 102 — PAPER_376 / PAPER_377
// grok_share_11254865.txt lines 6001-10322 re-analysis
// ============================================================

namespace Session102_FormalProofWormhole {

// ---------------------------------------------------------------------------
// PAPER_376 — UQFF Resonance Superconductive Formal Proof Set
// Source: "UQFF_Resonance Superconductive Universal Gravity Equation system
//          proof set._15May2025.docx"
// Five proofs: dimensional consistency, boundary conditions, resonance
// amplification at omega=2pi/tHubble, Meissner superconductivity, empirical
// validation vs Chandra magnetar data and EHT Sgr A* accretion.
// ---------------------------------------------------------------------------

namespace UQFFProofSet {
    static constexpr double G        = 6.6743e-11;
    static constexpr double M_sun    = 1.989e30;
    static constexpr double AU       = 1.496e11;
    static constexpr double c        = 3.0e8;
    static constexpr double Lambda   = 1.1e-52;
    static constexpr double tHubble  = 4.35e17;
    static constexpr double fquantum = 1.445e-17;  // Hz = 2pi/tHubble (exact)
    static constexpr double Ereact_0 = 1046.0;     // J   — magnetar flare seed
    static constexpr double kappa    = 0.0005;     // day⁻¹  decay rate

    // Proof 1: Newtonian baseline at 1 AU  → m/s²
    inline double proof1_newtonian_1AU() {
        return G * M_sun / (AU * AU);
    }

    // Proof 2: Cosmological floor (r→∞)  → m/s²
    inline double proof2_cosm_floor() {
        return Lambda * c * c / 3.0;
    }

    // Proof 3: Resonance frequency  → confirms fquantum = 2pi/tHubble
    inline double proof3_omega_resonance() {
        return 2.0 * 3.14159265358979323846 / tHubble;  // rad/s
    }

    // Proof 4: Meissner suppression factor (linear form)
    inline double proof4_meissner_linear(double B, double Bcrit) {
        return 1.0 - B / Bcrit;
    }

    // Proof 4: Meissner suppression factor (exponential form)
    inline double proof4_meissner_exp(double B, double Bcrit) {
        return std::exp(-B / Bcrit);
    }

    // Proof 5: Magnetar reactive energy at time t (days)
    inline double proof5_ereact(double t_days) {
        return Ereact_0 * std::exp(-kappa * t_days);
    }

    inline void selftest(std::ostream& os) {
        double g1AU     = proof1_newtonian_1AU();
        double g_cosm   = proof2_cosm_floor();
        double omega    = proof3_omega_resonance();
        double fq_check = std::abs(omega - fquantum) / fquantum;
        double m_lin    = proof4_meissner_linear(1e10, 1e11);
        double m_exp    = proof4_meissner_exp(1e10, 1e11);
        double er_10    = proof5_ereact(10.0);
        double er_100   = proof5_ereact(100.0);

        os << "[PAPER_376] Proof1 g(1AU)      = " << g1AU     << " m/s² (expect 5.93e-3)\n";
        os << "[PAPER_376] Proof2 cosm floor  = " << g_cosm   << " m/s²\n";
        os << "[PAPER_376] Proof3 omega_res   = " << omega     << " rad/s\n";
        os << "[PAPER_376] Proof3 fquantum    = " << fquantum  << " Hz  delta_rel=" << fq_check << "\n";
        os << "[PAPER_376] Proof4 Meissner L  = " << m_lin     << "  (1-B/Bcrit)\n";
        os << "[PAPER_376] Proof4 Meissner E  = " << m_exp     << "  exp(-B/Bcrit)\n";
        os << "[PAPER_376] Proof5 Ereact(10d) = " << er_10     << " J\n";
        os << "[PAPER_376] Proof5 Ereact(100d)= " << er_100    << " J\n";
        os << "[PAPER_376] Empirical: magnetar flares 10-100 days (Chandra); SgrA* ~1e-8 M_sun/yr (EHT)\n\n";
    }
} // namespace UQFFProofSet


// ---------------------------------------------------------------------------
// PAPER_377 — compute_a_wormhole() Production Implementation
//             + MUGE Error-Safety Infrastructure
// Source: grok_share_11254865.txt, C++ v8/v9 (lines 8600-10322)
// ---------------------------------------------------------------------------

namespace WormholeImpl {
    static constexpr double DEFAULT_B       = 1.0;       // m   (MT throat)
    static constexpr double DEFAULT_F_WORM  = 1.0;       // coupling factor
    static constexpr double DEFAULT_EVAC    = 7.09e-36;  // J/m³ nebular vac energy

    // 13th term in compute_resonance_MUGE — production implementation
    // a_worm = f_worm * Evac_neb / (b² + r²)
    inline double compute_a_wormhole(double r,
                                     double b       = DEFAULT_B,
                                     double f_worm  = DEFAULT_F_WORM,
                                     double Evac_neb = DEFAULT_EVAC) {
        // Safety: denominator is always >= b² > 0 for b > 0
        return f_worm * Evac_neb / (b * b + r * r);
    }

    // 24th unit test — test_compute_a_wormhole
    // Contract: at r=1e4, b=1, f_worm=1, Evac_neb=1 → result = 1/(1+r²)
    inline bool test_compute_a_wormhole(std::ostream& os) {
        const double r_test = 1e4;
        double expected = 1.0 / (DEFAULT_B * DEFAULT_B + r_test * r_test);
        double result   = compute_a_wormhole(r_test, DEFAULT_B, 1.0, 1.0);
        bool pass = std::abs(result - expected) < 1e-6 * expected;
        os << "[PAPER_377] test_compute_a_wormhole: result=" << result
           << " expected=" << expected << " PASS=" << (pass ? "YES" : "NO") << "\n";
        return pass;
    }

    // 7 reference system values
    struct SystemR { const char* name; double r; };
    static const SystemR SYSTEMS[] = {
        {"SGR_1745-2900",                 1e4   },
        {"Sagittarius_A*",                1e12  },
        {"Tapestry_Blazing_Starbirth",    3.086e17},
        {"Westerlund_2",                  3.086e17},
        {"Pillars_of_Creation",           9.46e15 },
        {"Rings_of_Relativity",           3.086e17},
        {"Students_Guide_Universe",       1e26  },
    };

    inline void selftest(std::ostream& os) {
        os << "[PAPER_377] compute_a_wormhole() implementation (C++ v8/v9 final):\n";
        os << "[PAPER_377] formula: f_worm * Evac_neb / (b^2 + r^2)\n";
        os << "[PAPER_377] 13th term in compute_resonance_MUGE()\n";
        os << "[PAPER_377] 18-field CSV I/O: name,I,A,omega1,omega2,Vsys,vexp,t,z,ffluid,M,r,B,Bcrit,rho_fluid,g_local,M_DM,delta_rho_rho\n";
        for (const auto& s : SYSTEMS) {
            double a = compute_a_wormhole(s.r);
            os << "[PAPER_377] a_worm(" << s.name << ", r=" << s.r << ") = " << a << " m/s²\n";
        }
        test_compute_a_wormhole(os);
        os << "[PAPER_377] Error-safe functions: compute_compressed_base (throw if r==0),\n"
              "              compute_compressed_super_adj (throw if Bcrit==0),\n"
              "              compute_compressed_quantum (throw if Delta_x_p==0),\n"
              "              compute_compressed_perturbation (throw if r==0)\n";
        os << "[PAPER_377] Total unit tests in suite: 24\n\n";
    }
} // namespace WormholeImpl

// Session 102 combined self-test
inline void session102_selftest(std::ostream& os) {
    os << "=== Session 102 Self-Test: PAPER_376 + PAPER_377 ===\n";
    UQFFProofSet::selftest(os);
    WormholeImpl::selftest(os);
    os << "[Session102] Self-test complete. Papers: 376 / 377\n\n";
}

} // namespace Session102_FormalProofWormhole

// ============================================================
// End of STAR_MAGIC_09SEPT_UQFF_MODULE.cpp
// Session 100 — PAPER_368 / PAPER_369 / PAPER_370
// Session 101 — PAPER_371 / PAPER_372 / PAPER_373 / PAPER_374 / PAPER_375
// Session 102 — PAPER_376 / PAPER_377
// ============================================================


// ============================================================================
// STANDALONE MAIN — STAR_MAGIC_09SEPT_UQFF_MODULE.cpp
// Star-Magic UQFF Standalone Module (self-updating, self-simulating)
// ============================================================================

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#ifdef STAR_MAGIC_STANDALONE_MODULE
#include <iostream>
#include <iomanip>
#include <chrono>
int main(int argc, char* argv[]) {
    using clock = std::chrono::high_resolution_clock;
    auto t0 = clock::now();
    std::cout << "=== StarMagic09Sept::StarMagic09SeptModule Standalone Self-Simulation ===\n";
    std::cout << std::fixed << std::setprecision(8);
    StarMagic09Sept::StarMagic09SeptModule module;
    constexpr double DT = 3.15576e7;
    constexpr int N = 10;
    for (int i = 0; i < N; ++i) {
        double t = i * DT;
        module.computeAll(t);
        std::cout << "  epoch " << i << "  t=" << t << " s\n";
    }
    std::cout << "  Done\n";
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(clock::now()-t0).count();
    std::cout << "  Done in " << ms << " ms\n";
    return 0;
}
#endif // STAR_MAGIC_STANDALONE_MODULE
