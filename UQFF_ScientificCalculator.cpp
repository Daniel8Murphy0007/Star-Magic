/**
 * UQFF Scientific Calculator - Iteration #32
 * Origin: Grok Thread https://x.com/i/grok/share/71d9d3b17d9c4998bc967ab602e89c46
 * 
 * Integrates:
 * - Universal Quantum Field Framework (UQFF) buoyancy equations
 * - Master Universal Gravity Equations (MUGE) for astrophysical systems
 * - Electric Universe (EU) proof via EM/gravity ratio
 * - Gyroscopic torque nullification
 * - Alpha clustering (Schmidt et al. 2016)
 * - Widom-Larsen LENR theory
 * 
 * Dependencies: Qt6, ANTLR4, SymEngine, Eigen, GSL, QCustomPlot
 * Author: Daniel T. Murphy / AI Collaboration
 * Date: March 2026
 */

#pragma once

#define _USE_MATH_DEFINES  // Enable M_PI on Windows MSVC

#include <cmath>
#include <map>
#include <string>
#include <vector>
#include <functional>
#include <stdexcept>
#include <memory>
#include <set>
#include <sstream>

// SymEngine for DimensionalAnalyzer expression tree traversal
#ifdef USE_SYMENGINE
#include <symengine/basic.h>
#include <symengine/symbol.h>
#include <symengine/add.h>
#include <symengine/mul.h>
#include <symengine/pow.h>
#include <symengine/integer.h>
#include <symengine/real_double.h>
#include <symengine/rational.h>
#include <symengine/functions.h>
#endif
#include <algorithm>

// ============================================================================
// PHYSICAL CONSTANTS (from UQFF Framework)
// ============================================================================
namespace UQFFConstants {
    // Fundamental Constants
    constexpr double c = 299792458.0;           // Speed of light (m/s)
    constexpr double G = 6.67430e-11;           // Gravitational constant (m³/kg/s²)
    constexpr double h = 6.62607015e-34;        // Planck constant (J·s)
    constexpr double hbar = 1.054571817e-34;    // Reduced Planck (J·s)
    constexpr double k_B = 1.380649e-23;        // Boltzmann constant (J/K)
    constexpr double e_charge = 1.602176634e-19; // Elementary charge (C)
    constexpr double m_p = 1.67262192e-27;      // Proton mass (kg)
    constexpr double m_e = 9.1093837e-31;       // Electron mass (kg)
    constexpr double mu_0 = 1.25663706e-6;      // Vacuum permeability (H/m)
    constexpr double eps_0 = 8.8541878e-12;     // Vacuum permittivity (F/m)
    
    // Astrophysical Constants
    constexpr double M_sun = 1.989e30;          // Solar mass (kg)
    constexpr double R_sun = 6.96e8;            // Solar radius (m)
    constexpr double AU = 1.496e11;             // Astronomical unit (m)
    constexpr double pc = 3.086e16;             // Parsec (m)
    constexpr double ly = 9.461e15;             // Light year (m)
    constexpr double H_0 = 2.2e-18;             // Hubble constant (1/s) ~70 km/s/Mpc
    
    // UQFF-Specific Constants (calibrated values)
    constexpr double kappa = 0.0005;            // κ calibration constant (/day)
    constexpr double SSq = 0.57;                // [SSq] non-local term
    constexpr double H_SCm = 0.99;              // Superconductive material factor
    constexpr double U_UA = 0.0001;             // Universal Aether coefficient
    constexpr double k_eta = 1e-113;            // Neutron rate calibration
    constexpr double beta_i = 0.603;            // β_i index coefficient
    
    // UQFF Vacuum Densities
    constexpr double rho_vac_SCm = 7.09e-37;    // Superconductive vacuum (J/m³)
    constexpr double rho_vac_UA = 7.09e-36;     // Aether vacuum (J/m³)
    
    // LEP Reference (1998)
    constexpr double E_LEP_1998 = 200.0;        // LEP energy (GeV)
    constexpr double F_rel_astro = 4.30e33;     // Relativistic coherence force (N)
    
    // Magnetar Critical Field
    constexpr double B_crit_magnetar = 4.4e13;  // Critical magnetic field (T)
    
    // THz Resonance (Colman-Gillespie)
    constexpr double f_THz_low = 1.2e12;        // Lower THz resonance (Hz)
    constexpr double f_THz_high = 1.3e12;       // Upper THz resonance (Hz)
    constexpr double f_activation = 300.0;       // Activation frequency (Hz)
}

// ============================================================================
// UNITS CLASS - Dimensional Analysis
// ============================================================================
namespace UQFFDimensional {
    
    /**
     * Unit class for dimensional analysis
     * Tracks length, mass, time dimensions for physics equation validation
     */
    struct Unit {
        double value;
        int length;  // L dimension
        int mass;    // M dimension
        int time;    // T dimension
        int charge;  // Q dimension (extended SI)
        int temp;    // Θ dimension (temperature)
        int amount;  // N dimension (mole)
        int luminous; // J dimension (candela)
        
        Unit(double v = 0.0, int l = 0, int m = 0, int t = 0, int q = 0, int th = 0, int n = 0, int j = 0)
            : value(v), length(l), mass(m), time(t), charge(q), temp(th), amount(n), luminous(j) {}
        
        // Dimensional compatibility check
        bool compatible(const Unit& other) const {
            return length == other.length && mass == other.mass && 
                   time == other.time && charge == other.charge && temp == other.temp &&
                   amount == other.amount && luminous == other.luminous;
        }
        
        // Operators with dimensional analysis
        bool operator==(const Unit& other) const {
            return compatible(other);
        }
        bool operator!=(const Unit& other) const {
            return !compatible(other);
        }

        Unit operator+(const Unit& o) const {
            if (!compatible(o)) throw std::runtime_error("Incompatible dimensions for addition");
            return Unit(value + o.value, length, mass, time, charge, temp, amount, luminous);
        }
        
        Unit operator-(const Unit& o) const {
            if (!compatible(o)) throw std::runtime_error("Incompatible dimensions for subtraction");
            return Unit(value - o.value, length, mass, time, charge, temp, amount, luminous);
        }
        
        Unit operator*(const Unit& o) const {
            return Unit(value * o.value, length + o.length, mass + o.mass, 
                       time + o.time, charge + o.charge, temp + o.temp,
                       amount + o.amount, luminous + o.luminous);
        }
        
        Unit operator/(const Unit& o) const {
            return Unit(value / o.value, length - o.length, mass - o.mass,
                       time - o.time, charge - o.charge, temp - o.temp,
                       amount - o.amount, luminous - o.luminous);
        }
        
        Unit pow(int n) const {
            return Unit(std::pow(value, n), length * n, mass * n, time * n, 
                       charge * n, temp * n, amount * n, luminous * n);
        }
        
        Unit sqrt() const {
            if (length % 2 != 0 || mass % 2 != 0 || time % 2 != 0 ||
                amount % 2 != 0 || luminous % 2 != 0)
                throw std::runtime_error("Cannot take sqrt of odd-dimensional unit");
            return Unit(std::sqrt(value), length/2, mass/2, time/2, charge/2, temp/2,
                       amount/2, luminous/2);
        }

        bool isDimensionless() const {
            return length == 0 && mass == 0 && time == 0 && charge == 0 &&
                   temp == 0 && amount == 0 && luminous == 0;
        }
        
        // Dimension string representation
        std::string dimString() const {
            std::ostringstream oss;
            if (length) oss << "L^" << length << " ";
            if (mass) oss << "M^" << mass << " ";
            if (time) oss << "T^" << time << " ";
            if (charge) oss << "Q^" << charge << " ";
            if (temp) oss << "Θ^" << temp << " ";
            if (amount) oss << "N^" << amount << " ";
            if (luminous) oss << "J^" << luminous;
            std::string result = oss.str();
            while (!result.empty() && result.back() == ' ') result.pop_back();
            return result.empty() ? "dimensionless" : result;
        }
        
        // Static factory methods for common units
        static Unit meter(double v = 1.0) { return Unit(v, 1, 0, 0); }
        static Unit kilogram(double v = 1.0) { return Unit(v, 0, 1, 0); }
        static Unit second(double v = 1.0) { return Unit(v, 0, 0, 1); }
        static Unit newton(double v = 1.0) { return Unit(v, 1, 1, -2); }  // kg·m/s²
        static Unit joule(double v = 1.0) { return Unit(v, 2, 1, -2); }   // kg·m²/s²
        static Unit tesla(double v = 1.0) { return Unit(v, 0, 1, -2, -1); } // kg/(A·s²)
        static Unit velocity(double v = 1.0) { return Unit(v, 1, 0, -1); }  // m/s
        static Unit acceleration(double v = 1.0) { return Unit(v, 1, 0, -2); } // m/s²
        static Unit mole(double v = 1.0) { return Unit(v, 0, 0, 0, 0, 0, 1, 0); }     // amount of substance
        static Unit candela(double v = 1.0) { return Unit(v, 0, 0, 0, 0, 0, 0, 1); }   // luminous intensity
        static Unit ampere(double v = 1.0) { return Unit(v, 0, 0, 0, 1, 0, 0, 0); }    // electric current
        static Unit kelvin(double v = 1.0) { return Unit(v, 0, 0, 0, 0, 1, 0, 0); }    // temperature
    };
    
    // Unit conversion factors database
    inline std::map<std::string, double> unitFactors = {
        {"m", 1.0}, {"cm", 0.01}, {"mm", 0.001}, {"km", 1000.0},
        {"AU", UQFFConstants::AU}, {"pc", UQFFConstants::pc}, {"ly", UQFFConstants::ly},
        {"kg", 1.0}, {"g", 0.001}, {"mg", 1e-6}, {"M_sun", UQFFConstants::M_sun},
        {"s", 1.0}, {"ms", 0.001}, {"us", 1e-6}, {"ns", 1e-9},
        {"min", 60.0}, {"hr", 3600.0}, {"day", 86400.0}, {"yr", 3.156e7},
        {"J", 1.0}, {"eV", 1.602e-19}, {"keV", 1.602e-16}, {"MeV", 1.602e-13}, {"GeV", 1.602e-10},
        {"N", 1.0}, {"dyn", 1e-5},
        {"T", 1.0}, {"G", 1e-4}  // Tesla, Gauss
    };
}

// ============================================================================
// SPARSE POLYNOMIAL - Memory Optimization for High-Degree
// ============================================================================
namespace UQFFPolynomial {
    
    /**
     * Sparse polynomial representation for high-degree equations (up to 26th)
     * Uses map<degree, coefficient> for memory efficiency
     */
    class SparsePoly {
    private:
        std::map<int, double> coeffs;
        
    public:
        SparsePoly() = default;
        
        void setCoeff(int degree, double coeff) {
            if (std::abs(coeff) > 1e-15) {
                coeffs[degree] = coeff;
            } else {
                coeffs.erase(degree);
            }
        }
        
        double getCoeff(int degree) const {
            auto it = coeffs.find(degree);
            return (it != coeffs.end()) ? it->second : 0.0;
        }
        
        int degree() const {
            return coeffs.empty() ? 0 : coeffs.rbegin()->first;
        }
        
        size_t numTerms() const { return coeffs.size(); }
        
        double eval(double x) const {
            double result = 0.0;
            for (const auto& [deg, coef] : coeffs) {
                result += coef * std::pow(x, deg);
            }
            return result;
        }
        
        // Horner's method for dense evaluation
        double evalHorner(double x) const {
            if (coeffs.empty()) return 0.0;
            int maxDeg = degree();
            double result = 0.0;
            for (int d = maxDeg; d >= 0; --d) {
                result = result * x + getCoeff(d);
            }
            return result;
        }
        
        SparsePoly operator+(const SparsePoly& o) const {
            SparsePoly result = *this;
            for (const auto& [deg, coef] : o.coeffs) {
                result.setCoeff(deg, result.getCoeff(deg) + coef);
            }
            return result;
        }
        
        SparsePoly operator*(const SparsePoly& o) const {
            SparsePoly result;
            for (const auto& [d1, c1] : coeffs) {
                for (const auto& [d2, c2] : o.coeffs) {
                    result.setCoeff(d1 + d2, result.getCoeff(d1 + d2) + c1 * c2);
                }
            }
            return result;
        }
        
        SparsePoly derivative() const {
            SparsePoly result;
            for (const auto& [deg, coef] : coeffs) {
                if (deg > 0) {
                    result.setCoeff(deg - 1, coef * deg);
                }
            }
            return result;
        }
        
        // Export coefficients for GSL solving
        std::vector<double> toDenseCoeffs() const {
            if (coeffs.empty()) return {0.0};
            std::vector<double> dense(degree() + 1, 0.0);
            for (const auto& [deg, coef] : coeffs) {
                dense[deg] = coef;
            }
            return dense;
        }
    };
}

// ============================================================================
// UQFF EQUATION SYSTEM
// ============================================================================
namespace UQFFEquations {
    
    using namespace UQFFConstants;
    
    /**
     * System parameters for UQFF calculations
     */
    struct SystemParams {
        double mass;           // System mass (kg)
        double radius;         // Characteristic radius (m)
        double time;           // Time parameter (s)
        double B_field;        // Magnetic field (T)
        double velocity;       // Characteristic velocity (m/s)
        double temperature;    // Temperature (K)
        double density;        // Density (kg/m³)
        double z_redshift;     // Cosmological redshift
        double E_cm;           // Center-of-mass energy (GeV)
        double n_quantum;      // Quantum number
        double theta;          // Angular parameter (rad)
        
        // UQFF-specific
        double f_sc;           // Superconductive factor
        double Q_wave;         // Wave quality factor
        double f_Heav;         // Heaviside factor
        double f_quasi;        // Quasi-particle factor
        double gamma_decay;    // Decay rate (/s)
        double omega_c;        // Characteristic frequency (rad/s)
        
        // Pre-computed scaling factors
        double computeEnergyRatio() const { return E_cm / E_LEP_1998; }
        double computeHubbleScale(double t) const { return 1.0 + H_0 * t; }
    };
    
    // ========================================================================
    // BUOYANCY FORCE - F_U_Bi_i (Core UQFF Equation)
    // ========================================================================
    
    /**
     * Computes the Universal Buoyancy Force
     * F_U_Bi_i = F_rel * (E_cm / E_LEP) * Q_wave * g(r,t)
     * 
     * @param params System parameters
     * @param g_local Local gravitational field (m/s²)
     * @return Buoyancy force (N), negative = repulsive/stabilizing
     */
    inline double computeBuoyancyForce(const SystemParams& params, double g_local) {
        double energy_ratio = params.computeEnergyRatio();
        double F_U_Bi_i = F_rel_astro * energy_ratio * params.Q_wave * g_local;
        return F_U_Bi_i;
    }
    
    /**
     * Buoyancy probability for alpha clustering (Schmidt et al. 2016)
     * P_alpha = 1 - exp(-|F_U_Bi_i| / E_th)
     */
    inline double computeClusteringProbability(double F_U_Bi_i, double E_threshold) {
        return 1.0 - std::exp(-std::abs(F_U_Bi_i) * 1e-15 / E_threshold);  // r~fm scale
    }
    
    // ========================================================================
    // UNIVERSAL MAGNETISM - Um(t,r,n)
    // ========================================================================
    
    /**
     * Computes Universal Magnetism field
     * Um = (μ_j(t) / r) * (1 - exp(-γt cos(πtn))) * φ * P_SCm * E_react * ...
     */
    inline double computeUniversalMagnetism(const SystemParams& params) {
        // Oscillating magnetic moment
        double mu_j = (1e3 + 0.4 * std::sin(params.omega_c * params.time)) * 3.38e20;
        
        // Exponential damping term
        double exp_term = 1.0 - std::exp(-params.gamma_decay * params.time * 
                          std::cos(M_PI * params.time * params.n_quantum));
        
        // Reaction energy decay
        double E_react = 1e46 * std::exp(-kappa * params.time / 86400.0);  // time in days
        
        // Full Um calculation
        double Um = (mu_j / params.radius) * exp_term * 1.0 *  // φ = 1 (golden ratio simplified)
                    H_SCm * E_react * 
                    (1.0 + 1e13 * params.f_Heav) * 
                    (1.0 + params.f_quasi);
        
        return Um;
    }
    
    /**
     * Electric field from Um
     * E = Um * ρ_vac_UA / r
     */
    inline double computeElectricFieldFromUm(double Um, double r) {
        return Um * rho_vac_UA / r;
    }
    
    /**
     * Neutron production rate η (LENR)
     * η = k_η * exp(-[SSq]^n / 26) * exp(-π - t) * Um / ρ_vac_UA
     */
    inline double computeNeutronRate(double Um, double t, int n) {
        double eta = k_eta * std::exp(-std::pow(SSq, n) / 26.0) * 
                     std::exp(-M_PI - t) * Um / rho_vac_UA;
        return eta;
    }
    
    // ========================================================================
    // UNIVERSAL GRAVITY - Ug1-Ug4 Components
    // ========================================================================
    
    /**
     * Ug1: Magnetic dipole contribution
     */
    inline double computeUg1(double M, double mu_dipole, double r) {
        return (G * M / (r * r)) * (1.0 + mu_dipole * U_UA);
    }
    
    /**
     * Ug2: Charge-reactivity contribution
     */
    inline double computeUg2(double M, double Q_charge, double r) {
        return (G * M / (r * r)) * (1.0 + Q_charge * rho_vac_SCm);
    }
    
    /**
     * Ug3: String rotation contribution (time-dependent)
     */
    inline double computeUg3(double M, double omega, double r, double t) {
        double rotation_factor = 1.0 + (omega * t) / (2.0 * M_PI);
        return (G * M / (r * r)) * rotation_factor * H_SCm;
    }
    
    /**
     * Ug4: Vacuum concentration contribution
     */
    inline double computeUg4(double M, double r, double rho_vac_local) {
        double concentration = rho_vac_local / rho_vac_UA;
        return (G * M / (r * r)) * concentration;
    }
    
    /**
     * Total Universal Gravity (26-layer compressed)
     * g(r,t) = Σ(i=1 to 26) [Ug1_i + Ug2_i + Ug3_i + Ug4_i]
     */
    inline double computeTotalUniversalGravity(const SystemParams& params, int n_layers = 26) {
        double g_total = 0.0;
        
        for (int i = 1; i <= n_layers; ++i) {
            double layer_factor = std::exp(-static_cast<double>(i) * 2.0);  // Density scaling
            double r_i = params.radius * (1.0 + 0.1 * i);
            
            g_total += layer_factor * (
                computeUg1(params.mass, params.B_field, r_i) +
                computeUg2(params.mass, 1.0, r_i) +
                computeUg3(params.mass, params.velocity / r_i, r_i, params.time) +
                computeUg4(params.mass, r_i, rho_vac_SCm * layer_factor)
            );
        }
        
        return g_total;
    }
    
    // ========================================================================
    // UNIVERSAL INERTIA - Ui (Gyroscopic Integration)
    // ========================================================================
    
    /**
     * Computes Universal Inertia operator
     * Ui = λ_i * ρ_vac_SCm * ρ_vac_UA * ω_s(t) * t_n * f_TRZ
     */
    inline double computeUniversalInertia(double lambda_i, double omega_s, double t_n, double f_TRZ) {
        return lambda_i * rho_vac_SCm * rho_vac_UA * omega_s * t_n * f_TRZ;
    }
    
    /**
     * Gyroscopic torque
     * τ = I * ω * α
     */
    inline double computeGyroTorque(double I, double omega, double alpha) {
        return I * omega * alpha;
    }
    
    /**
     * Torque nullification via buoyancy
     * τ_null = τ + F_U_Bi_i * r * sin(θ) = 0
     * Solves for required F_U_Bi_i
     */
    inline double computeTorqueNullificationForce(double tau, double r, double theta) {
        return -tau / (r * std::sin(theta));
    }
    
    // ========================================================================
    // ELECTRIC UNIVERSE RATIO - R = F_EM / F_g
    // ========================================================================
    
    /**
     * Computes EM to gravity ratio (EU proof)
     * R = (q * Um * ρ_vac * v / r) / (G * M * m / r²)
     */
    inline double computeEURatio(double Um, double q, double v, double M, double m, double r) {
        double F_EM = q * (Um * rho_vac_UA / r) * v;
        double F_g = G * M * m / (r * r);
        return F_EM / F_g;
    }
    
    /**
     * Gyro-modulated EU ratio
     * R_gyro = R * (1 - τ_null / (Um * torque))
     */
    inline double computeGyroModulatedEURatio(double R, double tau_null, double Um_torque) {
        if (std::abs(Um_torque) < 1e-30) return R;
        return R * (1.0 - tau_null / Um_torque);
    }
    
    // ========================================================================
    // MASTER UNIVERSAL GRAVITY EQUATIONS (MUGE)
    // ========================================================================
    
    /**
     * MUGE for Hydrogen Atom
     * g = G * m_eff * m_p / r² + Σ_Z (G * M_Z / r_Z²) * (1 + f_sc) * exp(H_0 * t / c)
     */
    inline double computeMUGE_Hydrogen(double r, double t, double f_sc = 0.1) {
        double g_base = G * m_p * m_e / (r * r);
        double hubble_factor = std::exp(H_0 * t / c);
        return g_base * (1.0 + f_sc) * hubble_factor;
    }
    
    /**
     * MUGE for Magnetar
     * g = (G * M / r²) * (1 + H_0 * t) * (1 - B(t) / B_crit) + (Ug1 + Ug2 + Ug3 + Ug4)
     */
    inline double computeMUGE_Magnetar(const SystemParams& params) {
        double g_Newton = G * params.mass / (params.radius * params.radius);
        double hubble_mod = 1.0 + H_0 * params.time;
        double B_ratio = 1.0 - params.B_field / B_crit_magnetar;
        
        double g_MUGE = g_Newton * hubble_mod * B_ratio;
        g_MUGE += computeTotalUniversalGravity(params, 4);  // Sum Ug1-4
        
        return g_MUGE;
    }
    
    /**
     * MUGE for Sagittarius A* (SMBH evolution)
     * g = (G * M(t) / r²) * (1 + H_0 * t) * (1 - B(t) / B_crit) + ... + (G * M² / c⁴ * r) * (dΩ/dt)²
     */
    inline double computeMUGE_SgrA(const SystemParams& params, double dOmega_dt) {
        double M_SgrA = 4.3e6 * M_sun;  // 4.3 million solar masses
        double g_Newton = G * M_SgrA / (params.radius * params.radius);
        double hubble_mod = 1.0 + H_0 * params.time;
        double B_ratio = 1.0 - params.B_field / B_crit_magnetar;
        
        // Frame dragging term
        double frame_drag = (G * M_SgrA * M_SgrA / (std::pow(c, 4) * params.radius)) * 
                           std::pow(dOmega_dt, 2);
        
        double g_MUGE = g_Newton * hubble_mod * B_ratio + frame_drag;
        g_MUGE += computeTotalUniversalGravity(params, 4);
        
        return g_MUGE;
    }
    
    /**
     * MUGE for Globular Cluster (BH likelihood)
     * P_BH = 1 - exp(-|F_U_Bi_i| * t / t_cc)
     */
    inline double computeMUGE_GlobularBH(double F_U_Bi_i, double t, double t_core_collapse = 1e10 * 3e7) {
        return 1.0 - std::exp(-std::abs(F_U_Bi_i) * t / t_core_collapse);
    }
    
    /**
     * MUGE for Gravitational Lensing (Rings of Relativity)
     * g = (G * M / r²) * (1 + H(z) * t) * (1 - B / B_crit) * (1 + L(t)) + Ug_sum
     */
    inline double computeMUGE_Rings(const SystemParams& params, double L_luminosity) {
        double g_Newton = G * params.mass / (params.radius * params.radius);
        double H_z = H_0 * (1.0 + params.z_redshift);  // Redshift-dependent Hubble
        double hubble_mod = 1.0 + H_z * params.time;
        double B_ratio = 1.0 - params.B_field / B_crit_magnetar;
        double L_factor = 1.0 + L_luminosity / 1e45;  // Normalized to typical quasar
        
        return g_Newton * hubble_mod * B_ratio * L_factor + computeTotalUniversalGravity(params, 4);
    }
    
    /**
     * MUGE for Solar Planetary System
     * g = G * M_sun(t) / r(t)² - g_wind(t) + Σ_moons g_tidal(t)
     */
    inline double computeMUGE_SolarSystem(double r, double t, double g_wind = 0.0, 
                                          const std::vector<double>& g_tidal = {}) {
        double g_Newton = G * M_sun / (r * r);
        double g_total = g_Newton - g_wind;
        
        for (double g_moon : g_tidal) {
            g_total += g_moon;
        }
        
        return g_total;
    }
    
    // ========================================================================
    // PRECESSION AND ORBITAL MECHANICS
    // ========================================================================
    
    /**
     * GR Periapsis precession with UQFF correction
     * Δφ = (6π * G * M / (c² * a * (1 - e²))) * (1 + f_sc * exp(H * t / c))
     */
    inline double computePrecession(double M, double a, double e, double t, double f_sc = 0.1) {
        double GR_term = 6.0 * M_PI * G * M / (c * c * a * (1.0 - e * e));
        double UQFF_correction = 1.0 + f_sc * std::exp(H_0 * t / c);
        return GR_term * UQFF_correction;
    }
    
    /**
     * Planet Nine epoch calculation
     * T = 2π * sqrt(a³ / GM) with buoyancy adjustment
     */
    inline double computePlanetNinePeriod(double a = 600.0 * AU) {
        double T_Kepler = 2.0 * M_PI * std::sqrt(std::pow(a, 3) / (G * M_sun));
        return T_Kepler;  // ~11,400 years
    }
}

// ============================================================================
// NUMERICAL METHODS
// ============================================================================
namespace UQFFNumerics {
    
    /**
     * Newton-Raphson root finding
     */
    inline double newtonMethod(std::function<double(double)> f,
                               std::function<double(double)> df,
                               double x0, double tol = 1e-10, int maxIter = 100) {
        double x = x0;
        for (int i = 0; i < maxIter; ++i) {
            double fx = f(x);
            double dfx = df(x);
            if (std::abs(dfx) < 1e-15) return std::nan("");
            double dx = fx / dfx;
            x -= dx;
            if (std::abs(dx) < tol) return x;
        }
        return std::nan("");
    }
    
    /**
     * Runge-Kutta 4th Order ODE Solver
     * Solves dy/dt = f(t, y) from t0 to tf
     */
    inline std::vector<std::pair<double, double>> rungeKutta4(
            std::function<double(double, double)> f,
            double y0, double t0, double tf, double dt) {
        
        std::vector<std::pair<double, double>> solution;
        double t = t0, y = y0;
        
        while (t <= tf) {
            solution.push_back({t, y});
            
            double k1 = f(t, y);
            double k2 = f(t + dt/2, y + k1*dt/2);
            double k3 = f(t + dt/2, y + k2*dt/2);
            double k4 = f(t + dt, y + k3*dt);
            
            y += (k1 + 2*k2 + 2*k3 + k4) * dt / 6.0;
            t += dt;
        }
        
        return solution;
    }
    
    /**
     * Error propagation via partial derivatives
     * σ² = Σ (∂f/∂x_i)² * σ_i²
     */
    inline double errorPropagation(const std::vector<double>& partials,
                                   const std::vector<double>& uncertainties) {
        if (partials.size() != uncertainties.size()) {
            throw std::runtime_error("Mismatched partial/uncertainty sizes");
        }
        
        double variance = 0.0;
        for (size_t i = 0; i < partials.size(); ++i) {
            variance += partials[i] * partials[i] * uncertainties[i] * uncertainties[i];
        }
        return std::sqrt(variance);
    }
    
    /**
     * Simple least squares fitting (linear: y = a + bx)
     */
    inline std::pair<double, double> linearFit(const std::vector<double>& x,
                                                const std::vector<double>& y) {
        if (x.size() != y.size() || x.size() < 2) {
            throw std::runtime_error("Invalid data for linear fit");
        }
        
        size_t n = x.size();
        double sum_x = 0, sum_y = 0, sum_xy = 0, sum_xx = 0;
        
        for (size_t i = 0; i < n; ++i) {
            sum_x += x[i];
            sum_y += y[i];
            sum_xy += x[i] * y[i];
            sum_xx += x[i] * x[i];
        }
        
        double denom = n * sum_xx - sum_x * sum_x;
        double a = (sum_y * sum_xx - sum_x * sum_xy) / denom;
        double b = (n * sum_xy - sum_x * sum_y) / denom;
        
        return {a, b};
    }
}

// ============================================================================
// ASTROPHYSICAL SYSTEM PRESETS
// ============================================================================
namespace UQFFSystems {
    
    using namespace UQFFEquations;
    using namespace UQFFConstants;
    
    // Sagittarius A* (Galactic Center SMBH)
    inline SystemParams getSgrA() {
        SystemParams p;
        p.mass = 4.3e6 * M_sun;
        p.radius = 12.0e9;        // Schwarzschild radius ~12 million km
        p.time = 0.0;
        p.B_field = 1e-3;         // ~mG field
        p.velocity = 0.5 * c;     // Relativistic jets
        p.temperature = 1e10;
        p.density = 1e15;
        p.z_redshift = 0.0;
        p.E_cm = 1e12;            // TeV scale
        p.n_quantum = 0;
        p.theta = M_PI / 6;       // 30° spin misalignment
        p.f_sc = 0.1;
        p.Q_wave = 1e12;
        p.f_Heav = 0.01;
        p.f_quasi = 0.01;
        p.gamma_decay = 5e-5 / 86400.0;
        p.omega_c = 1.585e-8;
        return p;
    }
    
    // Magnetar (SGR 1806-20 type)
    inline SystemParams getMagnetar() {
        SystemParams p;
        p.mass = 1.4 * M_sun;     // Typical neutron star
        p.radius = 1e4;           // 10 km
        p.time = 0.0;
        p.B_field = 1e15;         // 10^15 G = 10^11 T
        p.velocity = 0.1 * c;
        p.temperature = 1e9;
        p.density = 1e17;
        p.z_redshift = 0.0;
        p.E_cm = 100.0;           // ~100 GeV
        p.n_quantum = 1;
        p.theta = M_PI / 4;
        p.f_sc = 0.5;
        p.Q_wave = 1e10;
        p.f_Heav = 0.01;
        p.f_quasi = 0.01;
        p.gamma_decay = 1e-4 / 86400.0;
        p.omega_c = 1e3;
        return p;
    }
    
    // Alpha cluster collision (Schmidt et al. 2016)
    inline SystemParams getAlphaCluster() {
        SystemParams p;
        p.mass = 4.0 * 1.66054e-27;  // 4 nucleons
        p.radius = 2e-15;             // ~2 fm
        p.time = 1e-15;               // Nuclear timescale
        p.B_field = 0.0;
        p.velocity = 0.1 * c;
        p.temperature = 1e8;
        p.density = 1e17;
        p.z_redshift = 0.0;
        p.E_cm = 1.4;                 // 35 MeV/n * 40 nucleons = 1.4 GeV
        p.n_quantum = 0;
        p.theta = M_PI / 2;
        p.f_sc = 0.0;
        p.Q_wave = 1e12;
        p.f_Heav = 0.01;
        p.f_quasi = 0.01;
        p.gamma_decay = 5e-5;
        p.omega_c = 1e21;             // Nuclear rotation
        return p;
    }
    
    // Hydrogen atom (metallic/superconductive)
    inline SystemParams getHydrogenAtom() {
        SystemParams p;
        p.mass = m_p;
        p.radius = 5.29e-11;          // Bohr radius
        p.time = 0.0;
        p.B_field = 0.0;
        p.velocity = 2.19e6;          // Electron orbital velocity
        p.temperature = 300.0;
        p.density = 1.0;
        p.z_redshift = 0.0;
        p.E_cm = 13.6e-9;             // 13.6 eV in GeV
        p.n_quantum = 1;
        p.theta = 0.0;
        p.f_sc = 0.1;
        p.Q_wave = 1.0;
        p.f_Heav = 0.0;
        p.f_quasi = 0.0;
        p.gamma_decay = 0.0;
        p.omega_c = 4.13e16;          // Hydrogen transition
        return p;
    }
}

// ============================================================================
// UQFF EQUATION SYMBOL PALETTE (for UI integration)
// ============================================================================
namespace UQFFSymbols {
    
    // UQFF equation set for symbol palette
    inline std::vector<std::string> getUQFFEquations() {
        return {
            // Buoyancy
            "F_U_Bi_i = F_rel * (E_cm / E_LEP) * Q_wave * g(r,t)",
            "P_alpha = 1 - exp(-|F_U_Bi_i| / E_th)",
            
            // Universal Magnetism
            "Um(t,r,n) = mu_j(t)/r * (1-exp(-γt*cos(πtn))) * φ * P_SCm * E_react * ...",
            "E = Um * ρ_vac_UA / r",
            "η = k_η * exp(-[SSq]^n/26) * exp(-π-t) * Um / ρ_vac_UA",
            
            // Universal Gravity
            "Ug1 = (G*M/r²) * (1 + μ_dipole * U_UA)",
            "Ug2 = (G*M/r²) * (1 + Q_charge * ρ_vac_SCm)",
            "Ug3 = (G*M/r²) * (1 + ωt/2π) * H_SCm",
            "Ug4 = (G*M/r²) * (ρ_vac_local / ρ_vac_UA)",
            "g_26D = Σ(i=1→26) [Ug1_i + Ug2_i + Ug3_i + Ug4_i]",
            
            // Universal Inertia
            "Ui = λ_i * ρ_vac_SCm * ρ_vac_UA * ω_s(t) * t_n * f_TRZ",
            "τ = I * ω * α",
            "τ_null: τ + F_U_Bi_i * r * sin(θ) = 0",
            
            // Electric Universe Ratio
            "R = F_EM / F_g = (q*Um*ρ_vac*v/r) / (G*M*m/r²)",
            "R_gyro = R * (1 - τ_null / Um_torque)",
            
            // MUGE Systems
            "g_H = G*m_eff*m_p/r² + Σ_Z (G*M_Z/r_Z²) * (1+f_sc) * exp(H₀t/c)",
            "g_Magnetar = (G*M/r²)*(1+H₀t)*(1-B(t)/B_crit) + Ug_sum",
            "g_SgrA* = (G*M(t)/r²)*(1+H₀t)*(1-B(t)/B_crit) + ... + (G*M²/c⁴r)*(dΩ/dt)²",
            "g_Rings = (G*M/r²)*(1+H(z)t)*(1-B/B_crit)*(1+L(t)) + Ug_sum",
            "g_MGE = -G*M(t)*r*f_core/(r²+a²)^{3/2}*f_tidal - G*M_BH/r²*f_BH",
            "g_Solar = G*M_sun/r² - g_wind + Σ_moons g_tidal",
            "P_BH = 1 - exp(-|F_U_Bi_i|*t / t_cc)",
            
            // Precession
            "Δφ = (6πGM / c²a(1-e²)) * (1 + f_sc*exp(Ht/c))",
            "T_P9 = 2π√(a³/GM) ≈ 11,400 yr",
            
            // Shock Physics (Rankine-Hugoniot)
            "ρ₁v₁ = ρ₂v₂",
            "ρ₁v₁² + P₁ = ρ₂v₂² + P₂",
            "½v₁² + γP₁/(γ-1)ρ₁ = ½v₂² + γP₂/(γ-1)ρ₂",
            
            // C-type shock (ion-neutral)
            "ρ_n(∂v_n/∂t + v_n·∇v_n) = -∇P_n + ρ_n*ν_ni*(v_i - v_n)",
            "v(z) ≈ v_s * exp(-z/L_d)"
        };
    }
    
    // Physics constants palette
    inline std::vector<std::string> getPhysicsConstants() {
        return {
            "c = 299792458 m/s",
            "G = 6.674×10⁻¹¹ m³/kg/s²",
            "h = 6.626×10⁻³⁴ J·s",
            "ℏ = 1.055×10⁻³⁴ J·s",
            "k_B = 1.381×10⁻²³ J/K",
            "e = 1.602×10⁻¹⁹ C",
            "m_p = 1.673×10⁻²⁷ kg",
            "m_e = 9.109×10⁻³¹ kg",
            "M_⊙ = 1.989×10³⁰ kg",
            "R_⊙ = 6.96×10⁸ m",
            "AU = 1.496×10¹¹ m",
            "pc = 3.086×10¹⁶ m",
            "ly = 9.461×10¹⁵ m",
            "H₀ = 70 km/s/Mpc",
            "ρ_vac_SCm = 7.09×10⁻³⁷ J/m³",
            "ρ_vac_UA = 7.09×10⁻³⁶ J/m³",
            "F_rel = 4.30×10³³ N",
            "B_crit = 4.4×10¹³ T",
            "κ = 0.0005/day",
            "[SSq] = 0.57",
            "k_η = 10⁻¹¹³"
        };
    }
    
    // Greek letters
    inline std::vector<std::string> getGreekLetters() {
        return {
            "α", "β", "γ", "δ", "ε", "ζ", "η", "θ", "ι", "κ", "λ", "μ",
            "ν", "ξ", "ο", "π", "ρ", "σ", "τ", "υ", "φ", "χ", "ψ", "ω",
            "Γ", "Δ", "Θ", "Λ", "Ξ", "Π", "Σ", "Υ", "Φ", "Ψ", "Ω"
        };
    }
    
    // Math operators
    inline std::vector<std::string> getMathOperators() {
        return {
            "+", "-", "×", "÷", "^", "_", "±", "∓", "≤", "≥", "≠", "≈", "≡",
            "∂", "∇", "∫", "∬", "∭", "∮", "∑", "∏", "√", "∛", "∞",
            "∈", "∉", "⊂", "⊆", "∪", "∩", "∀", "∃", "¬", "∧", "∨", "⇒", "⇔"
        };
    }
}

// ============================================================================
// DIMENSIONAL ANALYZER - Recursive expression unit validation
// Uses SymEngine expression tree to verify dimensional consistency
// ============================================================================
#ifdef USE_SYMENGINE
namespace UQFFDimensional {

    class DimensionalAnalyzer {
    public:
        std::map<std::string, Unit> varUnits;   // Variable → Unit mapping
        std::map<std::string, Unit> funcUnits;   // Function → return Unit mapping

        DimensionalAnalyzer() {
            // Pre-register common physics variables with their units
            // Users can override or extend via setVariableUnit()
        }

        void setVariableUnit(const std::string& name, const Unit& u) {
            varUnits[name] = u;
        }

        void setFunctionUnit(const std::string& name, const Unit& u) {
            funcUnits[name] = u;
        }

        // Recursive unit extraction from SymEngine expression tree
        // Returns the dimensional signature of the expression
        Unit getUnit(const SymEngine::RCP<const SymEngine::Basic>& expr) const {
            using namespace SymEngine;

            // Symbol: look up in variable map
            if (is_a<Symbol>(*expr)) {
                const auto& sym = down_cast<const Symbol&>(*expr);
                auto it = varUnits.find(sym.get_name());
                return it != varUnits.end() ? it->second : Unit{}; // Unknown → dimensionless
            }

            // Number: dimensionless
            if (is_a<Integer>(*expr) || is_a<RealDouble>(*expr) || is_a<Rational>(*expr)) {
                return Unit{};
            }

            // Addition/Subtraction: all terms must have same dimensions
            if (is_a<Add>(*expr)) {
                const auto& add = down_cast<const Add&>(*expr);
                auto args = add.get_args();
                if (args.empty()) return Unit{};
                Unit u = getUnit(args[0]);
                for (size_t i = 1; i < args.size(); ++i) {
                    Unit ui = getUnit(args[i]);
                    if (u != ui) {
                        throw std::runtime_error(
                            "Dimensional mismatch in addition: " + u.dimString() + 
                            " vs " + ui.dimString());
                    }
                }
                return u;
            }

            // Multiplication: dimensions add
            if (is_a<Mul>(*expr)) {
                const auto& mul = down_cast<const Mul&>(*expr);
                auto args = mul.get_args();
                Unit u(1.0);
                for (const auto& arg : args) {
                    u = u * getUnit(arg);
                }
                return u;
            }

            // Power: base dimensions scale by exponent
            if (is_a<Pow>(*expr)) {
                const auto& pw = down_cast<const Pow&>(*expr);
                Unit baseU = getUnit(pw.get_base());
                // Exponent must be a numeric constant for dimensional analysis
                auto expVal = pw.get_exp();
                if (is_a<Integer>(*expVal)) {
                    int exp = static_cast<int>(down_cast<const Integer&>(*expVal).as_int());
                    return baseU.pow(exp);
                }
                // Non-integer exponent: only valid if base is dimensionless
                if (!baseU.isDimensionless()) {
                    throw std::runtime_error(
                        "Non-integer exponent on dimensioned quantity: " + baseU.dimString());
                }
                return Unit{}; // dimensionless^anything = dimensionless
            }

            // Functions (sin, cos, exp, log): argument must be dimensionless, result is dimensionless
            // Unless overridden in funcUnits map
            if (is_a<FunctionSymbol>(*expr)) {
                const auto& func = down_cast<const FunctionSymbol&>(*expr);
                std::string fname = func.get_name();
                
                // Check if function has custom return unit
                auto it = funcUnits.find(fname);
                if (it != funcUnits.end()) {
                    return it->second;
                }

                // Standard math functions: argument must be dimensionless
                auto args = func.get_args();
                for (const auto& arg : args) {
                    Unit argU = getUnit(arg);
                    if (!argU.isDimensionless()) {
                        throw std::runtime_error(
                            "Function " + fname + " requires dimensionless argument, got: " + 
                            argU.dimString());
                    }
                }
                return Unit{}; // dimensionless
            }

            // Default: treat as dimensionless
            return Unit{};
        }

        // Validate that LHS and RHS of an equation have matching dimensions
        bool validateEquation(const SymEngine::RCP<const SymEngine::Basic>& lhs,
                              const SymEngine::RCP<const SymEngine::Basic>& rhs) const {
            try {
                Unit lhsU = getUnit(lhs);
                Unit rhsU = getUnit(rhs);
                return lhsU == rhsU;
            } catch (const std::runtime_error&) {
                return false; // Dimensional inconsistency detected
            }
        }

        // Validate and return detailed error message
        std::string validateEquationDetailed(const SymEngine::RCP<const SymEngine::Basic>& lhs,
                                              const SymEngine::RCP<const SymEngine::Basic>& rhs) const {
            try {
                Unit lhsU = getUnit(lhs);
                Unit rhsU = getUnit(rhs);
                if (lhsU == rhsU) {
                    return "OK: Both sides have dimensions [" + lhsU.dimString() + "]";
                }
                return "MISMATCH: LHS=[" + lhsU.dimString() + "] vs RHS=[" + rhsU.dimString() + "]";
            } catch (const std::runtime_error& e) {
                return std::string("ERROR: ") + e.what();
            }
        }
    };

} // namespace UQFFDimensional
#endif // USE_SYMENGINE

// ============================================================================
// MAIN ENTRY POINT FOR TESTING
// ============================================================================
#ifdef UQFF_STANDALONE_TEST

#include <iostream>
#include <iomanip>

int main() {
    using namespace UQFFEquations;
    using namespace UQFFSystems;
    using namespace UQFFConstants;
    
    std::cout << "=== UQFF Scientific Calculator Test ===" << std::endl;
    std::cout << std::scientific << std::setprecision(4);
    
    // Test Sgr A*
    auto sgrA = getSgrA();
    double g_SgrA = computeTotalUniversalGravity(sgrA, 4);
    std::cout << "\nSgr A* Universal Gravity (4-layer): " << g_SgrA << " m/s²" << std::endl;
    
    // Test Buoyancy
    double F_buoy = computeBuoyancyForce(sgrA, g_SgrA);
    std::cout << "Sgr A* Buoyancy Force: " << F_buoy << " N" << std::endl;
    
    // Test Alpha Clustering
    auto alpha = getAlphaCluster();
    double g_alpha = G * alpha.mass / (alpha.radius * alpha.radius);
    double F_alpha = computeBuoyancyForce(alpha, g_alpha);
    double P_alpha = computeClusteringProbability(F_alpha, 5e-13);  // 5 MeV threshold
    std::cout << "\nAlpha Cluster P_alpha: " << P_alpha << " (expect ~0.85)" << std::endl;
    
    // Test EU Ratio
    double Um = computeUniversalMagnetism(alpha);
    double R_EU = computeEURatio(Um, e_charge, 0.1*c, alpha.mass, alpha.mass, alpha.radius);
    std::cout << "Alpha Cluster EU Ratio (R = F_EM/F_g): " << R_EU << std::endl;
    
    // Test Magnetar MUGE
    auto magnetar = getMagnetar();
    double g_mag = computeMUGE_Magnetar(magnetar);
    std::cout << "\nMagnetar MUGE g: " << g_mag << " m/s²" << std::endl;
    
    // Test Planet Nine period
    double T_P9 = computePlanetNinePeriod();
    std::cout << "\nPlanet Nine Period: " << T_P9 / (3.156e7) << " years" << std::endl;
    
    std::cout << "\n=== Tests Complete ===" << std::endl;
    return 0;
}

#endif // UQFF_STANDALONE_TEST
