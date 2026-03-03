/**
 * Equation Solver
 * Thread b6d9bc22 Priority 2 - Iteration #30-31 Complete Implementation
 * 
 * Multi-type equation solver integrating ANTLR4, SymEngine, GSL, and UQFF
 */

#pragma once

#include "antlr4_parser.h"
#include "symengine_wrapper.h"
#include "polynomial_solver.h"
#include "uqff_equations.h"
#include "dimensional_analysis.h"
#include <vector>
#include <map>
#include <string>

namespace CalculatorAdvanced {

/**
 * @brief Solution result for functional equations
 */
struct FunctionalSolution {
    std::map<std::string, double> solutions;  // Variable → value
    std::string simplified_form;              // Simplified equation
    std::vector<std::string> solution_steps;  // Step-by-step derivation
    bool success;
    std::string error_message;
};

/**
 * @brief Solution result for parametric curves
 */
struct ParametricSolution {
    std::string x_function;  // x(t)
    std::string y_function;  // y(t)
    std::string z_function;  // z(t) (optional)
    std::vector<double> t_values;
    std::vector<std::vector<double>> points;  // [x, y, z] coordinates
    bool success;
    std::string error_message;
};

/**
 * @brief Solution result for ODE systems
 */
struct ODESolution {
    std::vector<double> t_values;
    std::map<std::string, std::vector<double>> y_values;  // Variable → time series
    std::string method;  // "RK4", "Adams-Bashforth", etc.
    double max_error;
    bool success;
    std::string error_message;
};

/**
 * @brief Solution result for series expansions
 */
struct SeriesSolution {
    std::string series_expansion;  // Taylor/Fourier series
    std::vector<double> coefficients;
    int order;
    double convergence_radius;
    bool success;
    std::string error_message;
};

/**
 * @brief Multi-type equation solver
 */
class EquationSolver {
public:
    EquationSolver();
    
    /**
     * @brief Solve equation from string input
     * @param input User-entered equation (auto-detected type)
     * @return Appropriate solution type
     */
    FunctionalSolution solve(const std::string& input);
    
    /**
     * @brief Solve functional equation (e.g., "f(x,y) = x^2 + y^2 = 25")
     */
    FunctionalSolution solveFunctional(const ParsedEquation& eq);
    
    /**
     * @brief Solve parametric curves (e.g., "x(t) = cos(t), y(t) = sin(t)")
     */
    ParametricSolution solveParametric(
        const ParsedEquation& eq,
        double t_min = 0.0,
        double t_max = 10.0,
        int num_points = 1000
    );
    
    /**
     * @brief Solve ODE system (e.g., "dy/dx = -x*y")
     */
    ODESolution solveODE(
        const ParsedEquation& eq,
        const std::map<std::string, double>& initial_conditions,
        double t_min = 0.0,
        double t_max = 10.0,
        double step_size = 0.01
    );
    
    /**
     * @brief Compute series expansion
     */
    SeriesSolution solveSeries(
        const ParsedEquation& eq,
        const std::string& variable,
        double expansion_point = 0.0,
        int order = 10
    );
    
    /**
     * @brief Solve polynomial equation
     */
    std::vector<std::complex<double>> solvePolynomial(const ParsedEquation& eq);
    
    /**
     * @brief Solve UQFF equation with parameter substitution
     */
    FunctionalSolution solveUQFF(
        const std::string& uqff_name,
        const std::map<std::string, double>& parameters
    );
    
    /**
     * @brief Set symbolic solver backend
     */
    void setSymbolicSolver(SymbolicSolver* solver);
    
    /**
     * @brief Set polynomial solver backend
     */
    void setPolynomialSolver(PolynomialSolver* poly_solver);
    
    /**
     * @brief Set UQFF catalog
     */
    void setUQFFCatalog(UQFFEquationCatalog* catalog);
    
    /**
     * @brief Set dimensional checker
     */
    void setDimensionalSystem(DimensionalSystem* dims);
    
    /**
     * @brief Enable/disable dimensional checking
     */
    void setDimensionalCheckEnabled(bool enabled);
    
    /**
     * @brief Get last error message
     */
    std::string getLastError() const;
    
private:
    ANTLR4Parser parser_;
    SymbolicSolver* symbolic_solver_ = nullptr;
    PolynomialSolver* polynomial_solver_ = nullptr;
    UQFFEquationCatalog* uqff_catalog_ = nullptr;
    DimensionalSystem* dimensional_system_ = nullptr;
    bool dimensional_check_enabled_ = true;
    std::string last_error_;
    
    // Helper methods
    bool checkDimensionalConsistency(const ParsedEquation& eq);
    std::vector<std::string> generateSolutionSteps(const SymbolicExpression& expr);
};

/**
 * @brief Example usage:
 * 
 * EquationSolver solver;
 * 
 * // 1. Functional equation
 * auto result1 = solver.solve("x^2 + y^2 = 25");
 * // Solves for x and y, provides multiple solutions
 * // result1.solutions = {{"x", {-5, -4, ..., 4, 5}}, {"y", {0, 3, ..., 3, 0}}}
 * 
 * // 2. Parametric curve (helix)
 * auto result2 = solver.solve("x(t) = cos(t), y(t) = sin(t), z(t) = t");
 * // result2.points = [[1, 0, 0], [0.998, 0.062, 0.1], ...]
 * 
 * // 3. ODE system (harmonic oscillator)
 * auto parsed = solver.parser_.parse("dy/dx = -x");
 * auto result3 = solver.solveODE(parsed, {{"y", 1.0}}, 0, 10, 0.01);
 * // result3.y_values["y"] = [1.0, 0.999, 0.995, ...]
 * 
 * // 4. Series expansion
 * auto parsed_series = solver.parser_.parse("sin(x)");
 * auto result4 = solver.solveSeries(parsed_series, "x", 0, 10);
 * // result4.series_expansion = "x - x^3/6 + x^5/120 - x^7/5040 + ..."
 * // result4.coefficients = [0, 1, 0, -1/6, 0, 1/120, ...]
 * 
 * // 5. UQFF equation
 * auto result5 = solver.solveUQFF("F_U_Bi_i", {
 *     {"M", 1.989e30},
 *     {"r", 1e16},
 *     {"t", 3.156e7},
 *     {"theta", M_PI/2}
 * });
 * // result5.solutions["F_U_Bi_i"] = 4.30e33  // N
 * // result5.solution_steps = [
 * //   "Step 1: Substitute M = 1.989e30 kg",
 * //   "Step 2: Integrate Ug1 + Ug2 + ... over r",
 * //   "Step 3: F_U_Bi_i = 4.30e33 N"
 * // ]
 */

} // namespace CalculatorAdvanced
