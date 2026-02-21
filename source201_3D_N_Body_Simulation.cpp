/*
================================================================================
SOURCE201: 3D N-Body Simulation Module
================================================================================
Part of Star-Magic UQFF Framework
Created: 21 Feb 2026
Based on: 108_3D N-Body Simulation_cpp_14Jan2026.txt

Self-Expanding Framework 2.0-Enhanced Module
Features:
  - Leapfrog symplectic integrator (KDK variant)
  - Barnes-Hut octree O(N log N) optimization
  - Energy conservation monitoring
  - Dynamic term registration
  - Runtime parameter tuning
  - State persistence

Integration: SOURCE201 block in MAIN_1_CoAnQi.cpp
Physics Terms: NBodyGravitationalTerm, BarnesHutForceTerm, LeapfrogIntegratorTerm
================================================================================
*/

#include <iostream>
#include <vector>
#include <array>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <memory>
#include <string>
#include <map>
#include <functional>
#include <sstream>

namespace SOURCE201 {

// ============================================================================
// Physical Constants
// ============================================================================
constexpr double G_NEWTON = 6.67430e-11;      // m³ kg⁻¹ s⁻²
constexpr double SOLAR_MASS = 1.989e30;       // kg
constexpr double AU = 1.496e11;               // m
constexpr double PARSEC = 3.086e16;           // m
constexpr double LIGHT_YEAR = 9.461e15;       // m

// ============================================================================
// Particle Structure
// ============================================================================
struct Particle {
    double mass;
    std::array<double, 3> position;
    std::array<double, 3> velocity;
    std::array<double, 3> force;
    std::array<double, 3> velocity_half;  // For leapfrog: v_{n+1/2}
    int id;
    std::string name;
    
    Particle(double m, std::array<double, 3> pos, std::array<double, 3> vel, int particle_id = 0) 
        : mass(m), position(pos), velocity(vel), 
          force{0.0, 0.0, 0.0}, velocity_half{0.0, 0.0, 0.0},
          id(particle_id), name("Particle_" + std::to_string(particle_id)) {}
    
    double kinetic_energy() const {
        double v_sq = 0.0;
        for(int k = 0; k < 3; ++k) {
            v_sq += velocity[k] * velocity[k];
        }
        return 0.5 * mass * v_sq;
    }
    
    double distance_to(const Particle& other) const {
        double r_sq = 0.0;
        for(int k = 0; k < 3; ++k) {
            double dx = other.position[k] - position[k];
            r_sq += dx * dx;
        }
        return std::sqrt(r_sq);
    }
};

// ============================================================================
// Barnes-Hut Octree Node for O(N log N) optimization
// ============================================================================
struct OctreeNode {
    std::array<double, 3> center;            // Center of this octant
    double half_size;                         // Half the side length
    double total_mass;                        // Total mass in this node
    std::array<double, 3> center_of_mass;    // Center of mass
    std::unique_ptr<OctreeNode> children[8]; // 8 octants
    int particle_index;                       // -1 if internal node or empty
    bool is_leaf;
    
    OctreeNode(std::array<double, 3> c, double hs) 
        : center(c), half_size(hs), total_mass(0.0), 
          center_of_mass{0.0, 0.0, 0.0}, particle_index(-1), is_leaf(true) {
        for(int i = 0; i < 8; ++i) children[i] = nullptr;
    }
    
    int get_octant(const std::array<double, 3>& pos) const {
        int octant = 0;
        if(pos[0] >= center[0]) octant |= 1;
        if(pos[1] >= center[1]) octant |= 2;
        if(pos[2] >= center[2]) octant |= 4;
        return octant;
    }
    
    void insert(int idx, const std::vector<Particle>& particles) {
        if(is_leaf && particle_index == -1) {
            particle_index = idx;
            return;
        }
        
        if(is_leaf && particle_index >= 0) {
            is_leaf = false;
            int old_idx = particle_index;
            particle_index = -1;
            
            int oct_old = get_octant(particles[old_idx].position);
            double hs = half_size * 0.5;
            std::array<double, 3> child_center;
            child_center[0] = center[0] + ((oct_old & 1) ? hs : -hs);
            child_center[1] = center[1] + ((oct_old & 2) ? hs : -hs);
            child_center[2] = center[2] + ((oct_old & 4) ? hs : -hs);
            children[oct_old] = std::make_unique<OctreeNode>(child_center, hs);
            children[oct_old]->insert(old_idx, particles);
        }
        
        int oct_new = get_octant(particles[idx].position);
        if(!children[oct_new]) {
            double hs = half_size * 0.5;
            std::array<double, 3> child_center;
            child_center[0] = center[0] + ((oct_new & 1) ? hs : -hs);
            child_center[1] = center[1] + ((oct_new & 2) ? hs : -hs);
            child_center[2] = center[2] + ((oct_new & 4) ? hs : -hs);
            children[oct_new] = std::make_unique<OctreeNode>(child_center, hs);
        }
        children[oct_new]->insert(idx, particles);
    }
    
    void compute_mass_distribution(const std::vector<Particle>& particles) {
        if(is_leaf) {
            if(particle_index >= 0) {
                total_mass = particles[particle_index].mass;
                center_of_mass = particles[particle_index].position;
            }
            return;
        }
        
        total_mass = 0.0;
        center_of_mass = {0.0, 0.0, 0.0};
        
        for(int i = 0; i < 8; ++i) {
            if(children[i]) {
                children[i]->compute_mass_distribution(particles);
                double child_mass = children[i]->total_mass;
                total_mass += child_mass;
                for(int k = 0; k < 3; ++k) {
                    center_of_mass[k] += child_mass * children[i]->center_of_mass[k];
                }
            }
        }
        
        if(total_mass > 0.0) {
            for(int k = 0; k < 3; ++k) {
                center_of_mass[k] /= total_mass;
            }
        }
    }
    
    void compute_force_on_particle(Particle& p, double theta, double G, double softening) const {
        if(total_mass == 0.0) return;
        
        std::array<double, 3> r_vec;
        double r_squared = 0.0;
        for(int k = 0; k < 3; ++k) {
            r_vec[k] = center_of_mass[k] - p.position[k];
            r_squared += r_vec[k] * r_vec[k];
        }
        double r = std::sqrt(r_squared + softening * softening);
        
        if(r < softening) return;
        
        double s = 2.0 * half_size;
        if(is_leaf || (s / r < theta)) {
            double f_mag = G * p.mass * total_mass / (r_squared + softening * softening);
            for(int k = 0; k < 3; ++k) {
                p.force[k] += f_mag * r_vec[k] / r;
            }
        } else {
            for(int i = 0; i < 8; ++i) {
                if(children[i]) {
                    children[i]->compute_force_on_particle(p, theta, G, softening);
                }
            }
        }
    }
};

// ============================================================================
// Self-Expanding Framework 2.0: Dynamic Physics Term Base
// ============================================================================
class NBodyPhysicsTerm {
protected:
    std::string name_;
    std::string description_;
    bool enabled_;
    double learning_rate_;
    bool logging_enabled_;
    std::map<std::string, double> parameters_;
    
public:
    NBodyPhysicsTerm(const std::string& name, const std::string& desc)
        : name_(name), description_(desc), enabled_(true), 
          learning_rate_(0.001), logging_enabled_(false) {}
    
    virtual ~NBodyPhysicsTerm() = default;
    
    virtual double compute(const std::vector<Particle>& particles, 
                          const std::map<std::string, double>& params) const = 0;
    virtual bool validate() const = 0;
    
    // Self-Expanding Framework 2.0 methods
    void setParameter(const std::string& key, double value) { parameters_[key] = value; }
    double getParameter(const std::string& key) const {
        auto it = parameters_.find(key);
        return (it != parameters_.end()) ? it->second : 0.0;
    }
    
    void setEnabled(bool e) { enabled_ = e; }
    bool isEnabled() const { return enabled_; }
    
    void setLearningRate(double lr) { learning_rate_ = lr; }
    double getLearningRate() const { return learning_rate_; }
    
    void setLogging(bool log) { logging_enabled_ = log; }
    
    const std::string& getName() const { return name_; }
    const std::string& getDescription() const { return description_; }
    
    void log(const std::string& msg) const {
        if(logging_enabled_) {
            std::cout << "[" << name_ << "] " << msg << std::endl;
        }
    }
    
    virtual std::string exportState() const {
        std::ostringstream oss;
        oss << "Term: " << name_ << "\n";
        oss << "Enabled: " << (enabled_ ? "true" : "false") << "\n";
        oss << "LearningRate: " << learning_rate_ << "\n";
        for(const auto& p : parameters_) {
            oss << "Param_" << p.first << ": " << p.second << "\n";
        }
        return oss.str();
    }
};

// ============================================================================
// Concrete Physics Terms
// ============================================================================

// Newtonian Gravitational Force Term
class NBodyGravitationalTerm : public NBodyPhysicsTerm {
public:
    NBodyGravitationalTerm() 
        : NBodyPhysicsTerm("NBodyGravitationalTerm", 
            "Direct N-body gravitational force: F = G*m1*m2/r^2") {
        setParameter("G", G_NEWTON);
        setParameter("softening", 1e8);  // 100 km default
    }
    
    double compute(const std::vector<Particle>& particles,
                   const std::map<std::string, double>& params) const override {
        double G = getParameter("G");
        double softening = getParameter("softening");
        double total_potential = 0.0;
        
        for(size_t i = 0; i < particles.size(); ++i) {
            for(size_t j = i + 1; j < particles.size(); ++j) {
                double r_sq = 0.0;
                for(int k = 0; k < 3; ++k) {
                    double dx = particles[j].position[k] - particles[i].position[k];
                    r_sq += dx * dx;
                }
                double r = std::sqrt(r_sq + softening * softening);
                if(r > 0.0) {
                    total_potential -= G * particles[i].mass * particles[j].mass / r;
                }
            }
        }
        
        log("Computed gravitational potential: " + std::to_string(total_potential) + " J");
        return total_potential;
    }
    
    bool validate() const override {
        return getParameter("G") > 0 && getParameter("softening") >= 0;
    }
};

// Barnes-Hut Approximation Term
class BarnesHutForceTerm : public NBodyPhysicsTerm {
public:
    BarnesHutForceTerm()
        : NBodyPhysicsTerm("BarnesHutForceTerm",
            "Barnes-Hut tree approximation: O(N log N) force calculation") {
        setParameter("theta", 0.5);      // Opening angle
        setParameter("G", G_NEWTON);
        setParameter("softening", 1e8);
    }
    
    double compute(const std::vector<Particle>& particles,
                   const std::map<std::string, double>& params) const override {
        // Returns multipole acceptance ratio (quality metric)
        double theta = getParameter("theta");
        log("Barnes-Hut with theta = " + std::to_string(theta));
        return theta;  // Smaller = more accurate, larger = faster
    }
    
    bool validate() const override {
        double theta = getParameter("theta");
        return theta > 0.0 && theta < 2.0;  // Reasonable range
    }
};

// Leapfrog Integrator Term
class LeapfrogIntegratorTerm : public NBodyPhysicsTerm {
public:
    LeapfrogIntegratorTerm()
        : NBodyPhysicsTerm("LeapfrogIntegratorTerm",
            "Symplectic leapfrog (KDK): v_{n+1/2}=v_n+a*dt/2, x_{n+1}=x_n+v_{n+1/2}*dt") {
        setParameter("dt", 1e4);           // Default timestep
        setParameter("energy_tolerance", 0.01);  // 1% energy drift tolerance
    }
    
    double compute(const std::vector<Particle>& particles,
                   const std::map<std::string, double>& params) const override {
        // Compute total kinetic energy (validation metric)
        double T = 0.0;
        for(const auto& p : particles) {
            T += p.kinetic_energy();
        }
        log("Total kinetic energy: " + std::to_string(T) + " J");
        return T;
    }
    
    bool validate() const override {
        return getParameter("dt") > 0 && getParameter("energy_tolerance") > 0;
    }
};

// Energy Conservation Monitor Term
class EnergyConservationTerm : public NBodyPhysicsTerm {
private:
    mutable double initial_energy_;
    
public:
    EnergyConservationTerm()
        : NBodyPhysicsTerm("EnergyConservationTerm",
            "Energy conservation monitor: E = T + U, tracks ΔE/E₀"),
          initial_energy_(0.0) {
        setParameter("warning_threshold", 0.01);   // 1%
        setParameter("caution_threshold", 0.001);  // 0.1%
    }
    
    void setInitialEnergy(double E0) const { initial_energy_ = E0; }
    
    double compute(const std::vector<Particle>& particles,
                   const std::map<std::string, double>& params) const override {
        auto it = params.find("current_energy");
        if(it == params.end()) return 0.0;
        
        double current_energy = it->second;
        if(initial_energy_ == 0.0) {
            initial_energy_ = current_energy;
            return 0.0;
        }
        
        double relative_error = std::abs((current_energy - initial_energy_) / initial_energy_);
        
        if(relative_error > getParameter("warning_threshold")) {
            log("WARNING: Energy drift > " + 
                std::to_string(getParameter("warning_threshold") * 100) + "%");
        }
        
        return relative_error;
    }
    
    bool validate() const override {
        return getParameter("warning_threshold") > 0;
    }
};

// ============================================================================
// NBodySimulation Class - Main Simulation Engine
// ============================================================================
class NBodySimulation {
private:
    // Physical parameters
    double G_;
    double dt_;
    double softening_;
    double theta_;
    double expansion_factor_;
    
    // State
    double current_time_;
    double initial_energy_;
    bool leapfrog_initialized_;
    
    // Flags
    bool use_barnes_hut_;
    bool use_leapfrog_;
    bool enable_logging_;
    
    // Data
    std::vector<Particle> particles_;
    std::unique_ptr<OctreeNode> octree_root_;
    
    // Self-Expanding Framework 2.0
    std::map<std::string, double> dynamic_parameters_;
    std::vector<std::unique_ptr<NBodyPhysicsTerm>> registered_terms_;
    
public:
    // Constructor
    NBodySimulation(int num_particles = 5, double dt = 1e4)
        : G_(G_NEWTON), dt_(dt), softening_(1e8), theta_(0.5),
          expansion_factor_(1.01), current_time_(0.0), initial_energy_(0.0),
          leapfrog_initialized_(false), use_barnes_hut_(false),
          use_leapfrog_(true), enable_logging_(false)
    {
        // Register default physics terms
        registered_terms_.push_back(std::make_unique<NBodyGravitationalTerm>());
        registered_terms_.push_back(std::make_unique<BarnesHutForceTerm>());
        registered_terms_.push_back(std::make_unique<LeapfrogIntegratorTerm>());
        registered_terms_.push_back(std::make_unique<EnergyConservationTerm>());
        
        // Initialize particles
        initialize_particles(num_particles);
    }
    
    // ========================================================================
    // Initialization
    // ========================================================================
    void initialize_particles(int num_particles, unsigned int seed = 42) {
        std::srand(seed);
        particles_.clear();
        particles_.reserve(num_particles);
        
        for(int i = 0; i < num_particles; ++i) {
            double mass = random_range(1e26, 1e27);  // Small stellar masses
            std::array<double, 3> pos = {
                random_range(-1e11, 1e11),
                random_range(-1e11, 1e11),
                random_range(-1e11, 1e11)
            };
            std::array<double, 3> vel = {
                random_range(-1e3, 1e3),
                random_range(-1e3, 1e3),
                random_range(-1e3, 1e3)
            };
            particles_.emplace_back(mass, pos, vel, i);
        }
        
        leapfrog_initialized_ = false;
    }
    
    // ========================================================================
    // Force Computation
    // ========================================================================
    void compute_forces_direct() {
        for(auto& pi : particles_) {
            pi.force = {0.0, 0.0, 0.0};
            for(const auto& pj : particles_) {
                if(pi.id == pj.id) continue;
                
                std::array<double, 3> r_vec;
                double r_sq = 0.0;
                for(int k = 0; k < 3; ++k) {
                    r_vec[k] = pj.position[k] - pi.position[k];
                    r_sq += r_vec[k] * r_vec[k];
                }
                
                double r = std::sqrt(r_sq + softening_ * softening_);
                if(r > 0.0) {
                    double f_mag = G_ * pi.mass * pj.mass / (r_sq + softening_ * softening_);
                    for(int k = 0; k < 3; ++k) {
                        pi.force[k] += f_mag * r_vec[k] / r;
                    }
                }
            }
        }
    }
    
    void build_octree() {
        double min_c = -1e12, max_c = 1e12;
        for(const auto& p : particles_) {
            for(int k = 0; k < 3; ++k) {
                min_c = std::min(min_c, p.position[k]);
                max_c = std::max(max_c, p.position[k]);
            }
        }
        
        double hs = (max_c - min_c) * 0.6;
        double mid = (min_c + max_c) * 0.5;
        std::array<double, 3> center = {mid, mid, mid};
        
        octree_root_ = std::make_unique<OctreeNode>(center, hs);
        
        for(size_t i = 0; i < particles_.size(); ++i) {
            octree_root_->insert(static_cast<int>(i), particles_);
        }
        
        octree_root_->compute_mass_distribution(particles_);
    }
    
    void compute_forces_barnes_hut() {
        build_octree();
        for(auto& p : particles_) {
            p.force = {0.0, 0.0, 0.0};
            octree_root_->compute_force_on_particle(p, theta_, G_, softening_);
        }
    }
    
    void compute_forces() {
        if(use_barnes_hut_) {
            compute_forces_barnes_hut();
        } else {
            compute_forces_direct();
        }
    }
    
    // ========================================================================
    // Integration Methods
    // ========================================================================
    void update_euler(double dt) {
        for(auto& p : particles_) {
            for(int k = 0; k < 3; ++k) {
                double acc = p.force[k] / p.mass;
                p.velocity[k] += acc * dt;
                p.position[k] += p.velocity[k] * dt;
            }
        }
    }
    
    void initialize_leapfrog() {
        compute_forces();
        for(auto& p : particles_) {
            for(int k = 0; k < 3; ++k) {
                double acc = p.force[k] / p.mass;
                p.velocity_half[k] = p.velocity[k] + acc * dt_ * 0.5;
            }
        }
        leapfrog_initialized_ = true;
    }
    
    void update_leapfrog(double dt) {
        if(!leapfrog_initialized_) {
            initialize_leapfrog();
        }
        
        // Drift
        for(auto& p : particles_) {
            for(int k = 0; k < 3; ++k) {
                p.position[k] += p.velocity_half[k] * dt;
            }
        }
        
        // Force at new positions
        compute_forces();
        
        // Kick
        for(auto& p : particles_) {
            for(int k = 0; k < 3; ++k) {
                double acc = p.force[k] / p.mass;
                p.velocity[k] = p.velocity_half[k] + acc * dt * 0.5;
                p.velocity_half[k] += acc * dt;
            }
        }
    }
    
    void step(double dt) {
        if(use_leapfrog_) {
            update_leapfrog(dt);
        } else {
            compute_forces();
            update_euler(dt);
        }
        current_time_ += dt;
    }
    
    // ========================================================================
    // Energy Calculations
    // ========================================================================
    double compute_kinetic_energy() const {
        double T = 0.0;
        for(const auto& p : particles_) {
            T += p.kinetic_energy();
        }
        return T;
    }
    
    double compute_potential_energy() const {
        double U = 0.0;
        for(size_t i = 0; i < particles_.size(); ++i) {
            for(size_t j = i + 1; j < particles_.size(); ++j) {
                double r_sq = 0.0;
                for(int k = 0; k < 3; ++k) {
                    double dx = particles_[j].position[k] - particles_[i].position[k];
                    r_sq += dx * dx;
                }
                double r = std::sqrt(r_sq + softening_ * softening_);
                if(r > 0.0) {
                    U -= G_ * particles_[i].mass * particles_[j].mass / r;
                }
            }
        }
        return U;
    }
    
    double compute_total_energy() const {
        return compute_kinetic_energy() + compute_potential_energy();
    }
    
    void check_energy_conservation(int step) const {
        double E = compute_total_energy();
        double rel_err = std::abs((E - initial_energy_) / initial_energy_);
        
        std::cout << "Step " << step << ": E = " << E 
                  << " J, ΔE/E₀ = " << rel_err;
        
        if(rel_err > 0.01) {
            std::cout << " [WARNING: >1%]";
        } else if(rel_err > 0.001) {
            std::cout << " [Caution: >0.1%]";
        } else {
            std::cout << " [OK]";
        }
        std::cout << std::endl;
    }
    
    // ========================================================================
    // Self-Expanding Methods
    // ========================================================================
    void self_update() {
        dt_ *= 1.001;  // Adaptive timestep growth
    }
    
    void self_expand() {
        for(auto& p : particles_) {
            p.mass *= expansion_factor_;
            for(int k = 0; k < 3; ++k) {
                p.position[k] *= expansion_factor_;
            }
        }
    }
    
    // ========================================================================
    // Self-Expanding Framework 2.0 Interface
    // ========================================================================
    void registerDynamicTerm(std::unique_ptr<NBodyPhysicsTerm> term) {
        if(term->validate()) {
            registered_terms_.push_back(std::move(term));
        }
    }
    
    void setDynamicParameter(const std::string& key, double value) {
        dynamic_parameters_[key] = value;
    }
    
    double getDynamicParameter(const std::string& key) const {
        auto it = dynamic_parameters_.find(key);
        return (it != dynamic_parameters_.end()) ? it->second : 0.0;
    }
    
    std::string exportState() const {
        std::ostringstream oss;
        oss << "=== SOURCE201 NBodySimulation State ===\n";
        oss << "Time: " << current_time_ << " s\n";
        oss << "Particles: " << particles_.size() << "\n";
        oss << "dt: " << dt_ << " s\n";
        oss << "UseLeapfrog: " << (use_leapfrog_ ? "true" : "false") << "\n";
        oss << "UseBarnesHut: " << (use_barnes_hut_ ? "true" : "false") << "\n";
        oss << "Theta: " << theta_ << "\n";
        oss << "TotalEnergy: " << compute_total_energy() << " J\n";
        oss << "\n--- Registered Terms ---\n";
        for(const auto& term : registered_terms_) {
            oss << term->exportState() << "\n";
        }
        return oss.str();
    }
    
    void importState(const std::string& state) {
        // Parse and restore state (simplified)
        if(enable_logging_) {
            std::cout << "[SOURCE201] State import requested\n";
        }
    }
    
    void setEnableLogging(bool log) { 
        enable_logging_ = log;
        for(auto& term : registered_terms_) {
            term->setLogging(log);
        }
    }
    
    void setLearningRate(double lr) {
        for(auto& term : registered_terms_) {
            term->setLearningRate(lr);
        }
    }
    
    // ========================================================================
    // Main Simulation Loop
    // ========================================================================
    void simulate(int num_steps, bool output_files = true) {
        initial_energy_ = compute_total_energy();
        
        std::cout << "=== SOURCE201: 3D N-Body Simulation ===" << std::endl;
        std::cout << "Particles: " << particles_.size() << std::endl;
        std::cout << "Initial Energy: " << initial_energy_ << " J" << std::endl;
        std::cout << "Integrator: " << (use_leapfrog_ ? "LEAPFROG" : "EULER") << std::endl;
        std::cout << "Forces: " << (use_barnes_hut_ ? "BARNES-HUT O(N log N)" : "DIRECT O(N²)") << std::endl;
        std::cout << std::endl;
        
        for(int s = 0; s < num_steps; ++s) {
            step(dt_);
            
            if(output_files) {
                std::string fname = "nbody_step_" + std::to_string(s) + ".txt";
                std::ofstream out(fname);
                if(out) {
                    out << "# Step: " << s << ", Time: " << current_time_ << " s\n";
                    for(const auto& p : particles_) {
                        out << p.position[0] << " " << p.position[1] << " " << p.position[2] << "\n";
                    }
                }
            }
            
            if(s % 10 == 0) {
                check_energy_conservation(s);
            }
            
            self_update();
            if(s % 10 == 0) self_expand();
        }
        
        std::cout << "\n=== FINAL ENERGY REPORT ===" << std::endl;
        check_energy_conservation(num_steps);
    }
    
    // ========================================================================
    // Accessors
    // ========================================================================
    void set_use_leapfrog(bool use) { use_leapfrog_ = use; }
    void set_use_barnes_hut(bool use) { use_barnes_hut_ = use; }
    void set_theta(double t) { theta_ = t; }
    void set_softening(double s) { softening_ = s; }
    void set_dt(double dt) { dt_ = dt; }
    void set_G(double G) { G_ = G; }
    
    bool get_use_leapfrog() const { return use_leapfrog_; }
    bool get_use_barnes_hut() const { return use_barnes_hut_; }
    double get_theta() const { return theta_; }
    double get_current_time() const { return current_time_; }
    size_t get_num_particles() const { return particles_.size(); }
    const std::vector<Particle>& get_particles() const { return particles_; }
    
    // Add a particle manually
    void add_particle(double mass, std::array<double, 3> pos, std::array<double, 3> vel) {
        int id = static_cast<int>(particles_.size());
        particles_.emplace_back(mass, pos, vel, id);
    }
    
private:
    static double random_range(double min, double max) {
        return min + (std::rand() / static_cast<double>(RAND_MAX)) * (max - min);
    }
};

// ============================================================================
// UQFF Integration Functions (for MAIN_1_CoAnQi.cpp)
// ============================================================================

// Compute gravitational acceleration using direct summation
inline double compute_nbody_acceleration_SOURCE201(
    double M_total, double r, double softening = 1e8
) {
    return G_NEWTON * M_total / (r * r + softening * softening);
}

// Compute orbital period using Kepler's third law
inline double compute_orbital_period_SOURCE201(double a, double M_central) {
    return 2.0 * M_PI * std::sqrt(a * a * a / (G_NEWTON * M_central));
}

// Compute escape velocity
inline double compute_escape_velocity_SOURCE201(double M, double r) {
    return std::sqrt(2.0 * G_NEWTON * M / r);
}

// Compute virial temperature (for gas dynamics extension)
inline double compute_virial_temperature_SOURCE201(double M, double r) {
    constexpr double k_B = 1.380649e-23;  // Boltzmann constant J/K
    constexpr double m_H = 1.6735575e-27; // Hydrogen mass kg
    constexpr double mu = 0.6;            // Mean molecular weight
    return G_NEWTON * M * mu * m_H / (3.0 * k_B * r);
}

// Compute dynamical time
inline double compute_dynamical_time_SOURCE201(double rho) {
    return std::sqrt(3.0 * M_PI / (32.0 * G_NEWTON * rho));
}

// Compute Jeans mass (gravitational collapse threshold)
inline double compute_jeans_mass_SOURCE201(double T, double rho, double mu = 0.6) {
    constexpr double k_B = 1.380649e-23;
    constexpr double m_H = 1.6735575e-27;
    double c_s = std::sqrt(k_B * T / (mu * m_H));  // Sound speed
    return std::pow(M_PI / 6.0, 1.5) * std::pow(c_s, 3) / 
           std::sqrt(G_NEWTON * G_NEWTON * G_NEWTON * rho);
}

// ============================================================================
// Pre-defined Astrophysical Systems
// ============================================================================

// Solar System (inner planets)
inline NBodySimulation create_solar_system_inner_SOURCE201() {
    NBodySimulation sim(0, 86400.0);  // 1 day timestep
    
    // Sun at origin
    sim.add_particle(1.989e30, {0, 0, 0}, {0, 0, 0});
    
    // Mercury
    sim.add_particle(3.301e23, {5.79e10, 0, 0}, {0, 4.74e4, 0});
    
    // Venus
    sim.add_particle(4.867e24, {1.082e11, 0, 0}, {0, 3.50e4, 0});
    
    // Earth
    sim.add_particle(5.972e24, {1.496e11, 0, 0}, {0, 2.98e4, 0});
    
    // Mars
    sim.add_particle(6.417e23, {2.279e11, 0, 0}, {0, 2.41e4, 0});
    
    sim.set_softening(1e6);  // 1000 km softening
    return sim;
}

// Star cluster (Pleiades-like)
inline NBodySimulation create_star_cluster_SOURCE201(int num_stars = 100, double R = 3.0 * PARSEC) {
    NBodySimulation sim(0, 1e12);  // ~30,000 year timestep
    std::srand(12345);
    
    for(int i = 0; i < num_stars; ++i) {
        double r = R * std::pow(static_cast<double>(std::rand()) / RAND_MAX, 1.0/3.0);
        double theta = 2.0 * M_PI * std::rand() / RAND_MAX;
        double phi = std::acos(2.0 * std::rand() / static_cast<double>(RAND_MAX) - 1.0);
        
        std::array<double, 3> pos = {
            r * std::sin(phi) * std::cos(theta),
            r * std::sin(phi) * std::sin(theta),
            r * std::cos(phi)
        };
        
        double mass = SOLAR_MASS * (0.5 + 4.5 * std::rand() / static_cast<double>(RAND_MAX));
        
        // Velocity dispersion ~1 km/s
        double v_scale = 1e3;
        std::array<double, 3> vel = {
            v_scale * (2.0 * std::rand() / static_cast<double>(RAND_MAX) - 1.0),
            v_scale * (2.0 * std::rand() / static_cast<double>(RAND_MAX) - 1.0),
            v_scale * (2.0 * std::rand() / static_cast<double>(RAND_MAX) - 1.0)
        };
        
        sim.add_particle(mass, pos, vel);
    }
    
    sim.set_use_barnes_hut(true);  // Enable for large N
    sim.set_softening(1e14);       // ~0.003 pc softening
    return sim;
}

// Binary star system
inline NBodySimulation create_binary_star_SOURCE201(
    double M1 = SOLAR_MASS, double M2 = 0.8 * SOLAR_MASS,
    double separation = 1.0 * AU
) {
    NBodySimulation sim(0, 3600.0);  // 1 hour timestep
    
    double M_total = M1 + M2;
    double r1 = separation * M2 / M_total;
    double r2 = separation * M1 / M_total;
    
    double v_orb = std::sqrt(G_NEWTON * M_total / separation);
    double v1 = v_orb * M2 / M_total;
    double v2 = v_orb * M1 / M_total;
    
    sim.add_particle(M1, {-r1, 0, 0}, {0, -v1, 0});
    sim.add_particle(M2, {r2, 0, 0}, {0, v2, 0});
    
    sim.set_softening(1e8);
    return sim;
}

}  // namespace SOURCE201

// ============================================================================
// Standalone Execution (for testing)
// ============================================================================
#ifdef SOURCE201_STANDALONE

int main() {
    using namespace SOURCE201;
    
    std::cout << "=== SOURCE201: 3D N-Body Simulation Module ===" << std::endl;
    std::cout << "Star-Magic UQFF Framework" << std::endl;
    std::cout << "Self-Expanding Framework 2.0-Enhanced" << std::endl;
    std::cout << std::endl;
    
    // Create default simulation
    NBodySimulation sim(10, 1e4);
    sim.set_use_leapfrog(true);
    sim.set_use_barnes_hut(false);
    sim.set_theta(0.5);
    sim.setEnableLogging(true);
    
    // Run 50 steps
    sim.simulate(50, true);
    
    // Export state
    std::cout << "\n" << sim.exportState() << std::endl;
    
    // Test pre-defined systems
    std::cout << "\n=== Testing Solar System Inner ===" << std::endl;
    auto solar = create_solar_system_inner_SOURCE201();
    solar.set_use_leapfrog(true);
    solar.simulate(10, false);
    
    std::cout << "\nSOURCE201 module test complete." << std::endl;
    return 0;
}

#endif // SOURCE201_STANDALONE
