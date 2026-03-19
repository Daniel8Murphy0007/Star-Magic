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
