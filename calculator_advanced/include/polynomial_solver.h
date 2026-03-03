/**
 * Polynomial Solver (GSL 26th-Degree)
 * Thread b6d9bc22 Priority 2 - Iteration #30 Base Implementation
 * 
 * Solves polynomials up to degree 26 using GSL companion matrix method:
 * P(x) = a_n*x^n + a_(n-1)*x^(n-1) + ... + a_1*x + a_0 = 0
 */

#pragma once

#include <vector>
#include <complex>
#include <string>

namespace CalculatorAdvanced {

/**
 * @brief Polynomial coefficient representation
 */
struct Polynomial {
    std::vector<double> coefficients;  // [a_0, a_1, ..., a_n] (ascending powers)
    int degree;
    
    Polynomial(const std::vector<double>& coeffs);
    
    /**
     * @brief Evaluate polynomial at x
     */
    double evaluate(double x) const;
    std::complex<double> evaluate(std::complex<double> z) const;
    
    /**
     * @brief Get string representation
     */
    std::string toString() const;
};

/**
 * @brief Polynomial root finder using GSL
 */
class PolynomialSolver {
public:
    PolynomialSolver();
    ~PolynomialSolver();
    
    /**
     * @brief Solve polynomial equation P(x) = 0
     * @param poly Polynomial with coefficients
     * @return All roots (real and complex)
     */
    std::vector<std::complex<double>> solve(const Polynomial& poly);
    
    /**
     * @brief Find only real roots
     * @param poly Polynomial
     * @param tolerance Imaginary part threshold (default 1e-10)
     * @return Real roots only
     */
    std::vector<double> findRealRoots(const Polynomial& poly, double tolerance = 1e-10);
    
    /**
     * @brief Find roots in specific range [min, max]
     */
    std::vector<double> findRootsInRange(const Polynomial& poly, double min, double max);
    
    /**
     * @brief Verify root accuracy
     * @param poly Polynomial
     * @param root Candidate root
     * @return |P(root)| residual
     */
    double verifyRoot(const Polynomial& poly, std::complex<double> root) const;
    
    /**
     * @brief Get condition number (stability indicator)
     */
    double getConditionNumber() const;
    
private:
    // GSL workspace (gsl_poly_complex_workspace)
    void* workspace_ = nullptr;
    double condition_number_ = 1.0;
};

/**
 * @brief Example usage:
 * 
 * // Quadratic: x^2 - 5x + 6 = 0 (roots: 2, 3)
 * Polynomial quad({6, -5, 1});  // [a_0, a_1, a_2]
 * PolynomialSolver solver;
 * auto roots = solver.solve(quad);
 * // roots[0] = 2.0 + 0i, roots[1] = 3.0 + 0i
 * 
 * // 26th-degree UQFF vacuum energy polynomial
 * // P(E) = Σ(i=0 to 26) [c_i * (E/E_0)^i]
 * std::vector<double> coeffs_26D(27);  // 27 coefficients for degree 26
 * for (int i = 0; i <= 26; ++i) {
 *     coeffs_26D[i] = pow(10, -i*2);  // UQFF quantum level scaling
 * }
 * Polynomial uqff_poly(coeffs_26D);
 * auto uqff_roots = solver.solve(uqff_poly);
 * 
 * // Find only real physical roots
 * auto real_roots = solver.findRealRoots(uqff_poly);
 */

} // namespace CalculatorAdvanced
