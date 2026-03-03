#include "../include/polynomial_solver.h"
#include <gsl/gsl_poly.h>
#include <cmath>
#include <stdexcept>

Polynomial::Polynomial(const std::vector<double>& coeffs)
    : coefficients(coeffs), degree(coeffs.size() - 1)
{
}

PolynomialSolver::PolynomialSolver() {}

std::vector<std::complex<double>> PolynomialSolver::solve(const Polynomial& poly) {
    std::vector<std::complex<double>> roots;
    
    if (poly.degree > 26) {
        throw std::runtime_error("Polynomial degree exceeds maximum (26)");
    }
    
    // Use GSL for high-degree polynomials
    gsl_poly_complex_workspace* w = gsl_poly_complex_workspace_alloc(poly.degree + 1);
    std::vector<double> z(2 * poly.degree);
    
    gsl_poly_complex_solve(poly.coefficients.data(), poly.degree + 1, w, z.data());
    
    for (size_t i = 0; i < poly.degree; i++) {
        roots.push_back(std::complex<double>(z[2*i], z[2*i+1]));
    }
    
    gsl_poly_complex_workspace_free(w);
    return roots;
}

std::vector<double> PolynomialSolver::solveReal(const Polynomial& poly) {
    auto complexRoots = solve(poly);
    std::vector<double> realRoots;
    
    for (const auto& root : complexRoots) {
        if (std::abs(root.imag()) < 1e-10) {
            realRoots.push_back(root.real());
        }
    }
    
    return realRoots;
}
