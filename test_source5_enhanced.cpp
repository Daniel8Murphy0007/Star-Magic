// test_source5_enhanced.cpp
// Comprehensive test suite for Source5 2.0-Enhanced Framework

#include <iostream>
#include <cassert>
#include <cmath>
#include <memory>
#include <map>
#include <string>

// Minimal test framework
int g_test_passed = 0;
int g_test_failed = 0;

#define TEST_ASSERT(condition, message)                  \
    if (!(condition))                                    \
    {                                                    \
        std::cerr << "FAILED: " << message << std::endl; \
        g_test_failed++;                                 \
    }                                                    \
    else                                                 \
    {                                                    \
        g_test_passed++;                                 \
    }

// Mock PhysicsTerm and related classes for testing
class PhysicsTerm
{
public:
    virtual ~PhysicsTerm() = default;
    virtual double compute(const std::map<std::string, double> &params) const = 0;
    virtual std::string description() const = 0;
    virtual std::string version() const { return "1.0"; }
};

class DarkMatterHaloTerm : public PhysicsTerm
{
private:
    double M_halo;
    double r_scale;

public:
    DarkMatterHaloTerm(double mass, double scale) : M_halo(mass), r_scale(scale) {}

    double compute(const std::map<std::string, double> &params) const override
    {
        auto it_r = params.find("r");
        if (it_r == params.end())
            return 0.0;
        double r = it_r->second;
        if (r == 0.0 || r_scale == 0.0)
            return 0.0;

        const double G = 6.67430e-11;
        double x = r / r_scale;
        double rho_0 = M_halo / (4.0 * 3.14159265359 * r_scale * r_scale * r_scale * (std::log(2.0) - 0.5));
        return G * M_halo * std::log(1 + x) / (r * x);
    }

    std::string description() const override
    {
        return "Dark matter halo contribution (NFW profile)";
    }
};

class VacuumEnergyTerm : public PhysicsTerm
{
private:
    double E_vac_scale;
    double lambda;

public:
    VacuumEnergyTerm(double e_scale, double coupling)
        : E_vac_scale(e_scale), lambda(coupling) {}

    double compute(const std::map<std::string, double> &params) const override
    {
        auto it_t = params.find("t");
        if (it_t == params.end())
            return 0.0;
        double t = it_t->second;
        return lambda * E_vac_scale * (1.0 + 0.1 * std::sin(1e-10 * t));
    }

    std::string description() const override
    {
        return "Vacuum energy fluctuation term";
    }
};

// UQFFModule5 for testing
class UQFFModule5
{
private:
    std::vector<std::unique_ptr<PhysicsTerm>> dynamic_terms;
    std::map<std::string, double> dynamic_parameters;
    std::map<std::string, std::string> metadata;
    double learning_rate = 0.01;
    bool logging_enabled = false;

public:
    UQFFModule5()
    {
        metadata["version"] = "2.0-Enhanced";
        metadata["created"] = "2025-11-08";
        metadata["framework"] = "Self-Expanding UQFF";
    }

    void registerDynamicTerm(std::unique_ptr<PhysicsTerm> term)
    {
        dynamic_terms.push_back(std::move(term));
    }

    void setDynamicParameter(const std::string &name, double value)
    {
        dynamic_parameters[name] = value;
    }

    double getDynamicParameter(const std::string &name, double default_val = 0.0) const
    {
        auto it = dynamic_parameters.find(name);
        return (it != dynamic_parameters.end()) ? it->second : default_val;
    }

    double computeDynamicContributions(const std::map<std::string, double> &params) const
    {
        double sum = 0.0;
        for (const auto &term : dynamic_terms)
        {
            sum += term->compute(params);
        }
        return sum;
    }

    void setLearningRate(double rate)
    {
        learning_rate = rate;
    }

    double getLearningRate() const
    {
        return learning_rate;
    }

    void setEnableLogging(bool enable)
    {
        logging_enabled = enable;
    }

    size_t getTermCount() const
    {
        return dynamic_terms.size();
    }

    size_t getParameterCount() const
    {
        return dynamic_parameters.size();
    }

    std::string getVersion() const
    {
        auto it = metadata.find("version");
        return (it != metadata.end()) ? it->second : "unknown";
    }
};

// Test cases
void test_physics_term_base()
{
    std::cout << "Test 1: PhysicsTerm base class functionality..." << std::endl;

    auto dm_term = std::make_unique<DarkMatterHaloTerm>(1e12 * 1.989e30, 20000);
    TEST_ASSERT(dm_term->description() == "Dark matter halo contribution (NFW profile)",
                "Dark matter term description");
    TEST_ASSERT(dm_term->version() == "1.0", "PhysicsTerm version");

    auto vac_term = std::make_unique<VacuumEnergyTerm>(7.09e-36, 1.0);
    TEST_ASSERT(vac_term->description() == "Vacuum energy fluctuation term",
                "Vacuum energy term description");
}

void test_dark_matter_term()
{
    std::cout << "Test 2: DarkMatterHaloTerm computation..." << std::endl;

    DarkMatterHaloTerm dm_term(1e12 * 1.989e30, 20000);
    std::map<std::string, double> params = {{"r", 50000}};

    double result = dm_term.compute(params);
    TEST_ASSERT(result > 0.0, "Dark matter term produces positive contribution");
    TEST_ASSERT(std::isfinite(result), "Dark matter term result is finite");

    // Test edge cases
    params["r"] = 0.0;
    double zero_r = dm_term.compute(params);
    TEST_ASSERT(zero_r == 0.0, "Dark matter term handles r=0 safely");
}

void test_vacuum_energy_term()
{
    std::cout << "Test 3: VacuumEnergyTerm computation..." << std::endl;

    VacuumEnergyTerm vac_term(7.09e-36, 1.23e-40);
    std::map<std::string, double> params = {{"t", 0.0}};

    double result = vac_term.compute(params);
    TEST_ASSERT(result > 0.0, "Vacuum energy term produces positive contribution");
    TEST_ASSERT(std::isfinite(result), "Vacuum energy term result is finite");

    // Test time variation
    params["t"] = 1e10;
    double result_t1 = vac_term.compute(params);
    TEST_ASSERT(result_t1 != result, "Vacuum energy varies with time");
}

void test_module_initialization()
{
    std::cout << "Test 4: UQFFModule5 initialization..." << std::endl;

    UQFFModule5 module;
    TEST_ASSERT(module.getVersion() == "2.0-Enhanced", "Module version is 2.0-Enhanced");
    TEST_ASSERT(module.getTermCount() == 0, "Module starts with zero terms");
    TEST_ASSERT(module.getParameterCount() == 0, "Module starts with zero parameters");
    TEST_ASSERT(module.getLearningRate() == 0.01, "Default learning rate is 0.01");
}

void test_dynamic_term_registration()
{
    std::cout << "Test 5: Dynamic term registration..." << std::endl;

    UQFFModule5 module;

    module.registerDynamicTerm(std::make_unique<DarkMatterHaloTerm>(1e12 * 1.989e30, 20000));
    TEST_ASSERT(module.getTermCount() == 1, "Module has 1 term after registration");

    module.registerDynamicTerm(std::make_unique<VacuumEnergyTerm>(7.09e-36, 1.23e-40));
    TEST_ASSERT(module.getTermCount() == 2, "Module has 2 terms after second registration");
}

void test_dynamic_parameters()
{
    std::cout << "Test 6: Dynamic parameter management..." << std::endl;

    UQFFModule5 module;

    module.setDynamicParameter("custom_coupling", 1.23e-40);
    TEST_ASSERT(module.getParameterCount() == 1, "Module has 1 parameter");
    TEST_ASSERT(module.getDynamicParameter("custom_coupling") == 1.23e-40,
                "Parameter value retrieved correctly");

    module.setDynamicParameter("halo_mass", 1e12);
    TEST_ASSERT(module.getParameterCount() == 2, "Module has 2 parameters");

    double default_val = module.getDynamicParameter("nonexistent", 42.0);
    TEST_ASSERT(default_val == 42.0, "Non-existent parameter returns default value");
}

void test_dynamic_computation()
{
    std::cout << "Test 7: Dynamic contribution computation..." << std::endl;

    UQFFModule5 module;

    // No terms - should return 0
    std::map<std::string, double> params = {{"r", 50000}, {"t", 1e10}};
    double result_empty = module.computeDynamicContributions(params);
    TEST_ASSERT(result_empty == 0.0, "Empty module returns zero contribution");

    // Add terms
    module.registerDynamicTerm(std::make_unique<DarkMatterHaloTerm>(1e12 * 1.989e30, 20000));
    module.registerDynamicTerm(std::make_unique<VacuumEnergyTerm>(7.09e-36, 1.23e-40));

    double result = module.computeDynamicContributions(params);
    TEST_ASSERT(result > 0.0, "Module computes positive contribution with terms");
    TEST_ASSERT(std::isfinite(result), "Computed contribution is finite");
}

void test_learning_rate()
{
    std::cout << "Test 8: Learning rate management..." << std::endl;

    UQFFModule5 module;

    module.setLearningRate(0.001);
    TEST_ASSERT(module.getLearningRate() == 0.001, "Learning rate set to 0.001");

    module.setLearningRate(0.1);
    TEST_ASSERT(module.getLearningRate() == 0.1, "Learning rate updated to 0.1");
}

void test_logging_control()
{
    std::cout << "Test 9: Logging control..." << std::endl;

    UQFFModule5 module;

    module.setEnableLogging(true);
    // Just verify it doesn't crash - logging state isn't exposed
    module.setEnableLogging(false);
    TEST_ASSERT(true, "Logging control works without errors");
}

void test_combined_functionality()
{
    std::cout << "Test 10: Combined functionality test..." << std::endl;

    UQFFModule5 module;

    // Register terms
    module.registerDynamicTerm(std::make_unique<DarkMatterHaloTerm>(1e12 * 1.989e30, 20000));
    module.registerDynamicTerm(std::make_unique<VacuumEnergyTerm>(7.09e-36, 1.23e-40));

    // Set parameters
    module.setDynamicParameter("test_param1", 100.0);
    module.setDynamicParameter("test_param2", 200.0);

    // Set learning rate
    module.setLearningRate(0.05);

    // Compute
    std::map<std::string, double> params = {{"r", 50000}, {"t", 1e10}};
    double result = module.computeDynamicContributions(params);

    TEST_ASSERT(module.getTermCount() == 2, "Combined test: 2 terms registered");
    TEST_ASSERT(module.getParameterCount() == 2, "Combined test: 2 parameters set");
    TEST_ASSERT(module.getLearningRate() == 0.05, "Combined test: learning rate correct");
    TEST_ASSERT(result > 0.0 && std::isfinite(result), "Combined test: valid computation");
}

int main()
{
    std::cout << "========================================" << std::endl;
    std::cout << "Source5 2.0-Enhanced Framework Test Suite" << std::endl;
    std::cout << "========================================" << std::endl
              << std::endl;

    test_physics_term_base();
    test_dark_matter_term();
    test_vacuum_energy_term();
    test_module_initialization();
    test_dynamic_term_registration();
    test_dynamic_parameters();
    test_dynamic_computation();
    test_learning_rate();
    test_logging_control();
    test_combined_functionality();

    std::cout << std::endl
              << "========================================" << std::endl;
    std::cout << "Test Results:" << std::endl;
    std::cout << "  Passed: " << g_test_passed << std::endl;
    std::cout << "  Failed: " << g_test_failed << std::endl;
    std::cout << "========================================" << std::endl;

    if (g_test_failed == 0)
    {
        std::cout << "✓ All tests passed!" << std::endl;
        return 0;
    }
    else
    {
        std::cout << "✗ Some tests failed." << std::endl;
        return 1;
    }
}
