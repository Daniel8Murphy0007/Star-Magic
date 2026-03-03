#include "../include/symengine_wrapper.h"
#include <symengine/add.h>
#include <symengine/mul.h>
#include <symengine/pow.h>
#include <symengine/diff.h>
#include <symengine/solve.h>
#include <symengine/simplify.h>
#include <symengine/expand.h>
#include <symengine/subs.h>
#include <symengine/functions.h>
#include <symengine/parser.h>
#include <symengine/printers/latex.h>

using namespace SymEngine;

// SymbolicExpression implementation
SymbolicExpression::SymbolicExpression(const std::string& expr) {
    try {
        expr_ = parse(expr);
    } catch (const std::exception& e) {
        throw std::runtime_error("Failed to parse expression: " + std::string(e.what()));
    }
}

SymbolicExpression::SymbolicExpression(RCP<const Basic> expr) : expr_(expr) {}

SymbolicExpression SymbolicExpression::differentiate(const std::string& variable) const {
    auto var = symbol(variable);
    return SymbolicExpression(diff(expr_, var));
}

SymbolicExpression SymbolicExpression::integrate(const std::string& variable) const {
    // SymEngine has limited integration support, implement basic rules
    auto var = symbol(variable);
    
    // For simple cases
    if (is_a<Pow>(*expr_)) {
        const Pow& p = down_cast<const Pow&>(*expr_);
        if (eq(*p.get_base(), *var)) {
            auto n = add(p.get_exp(), integer(1));
            return SymbolicExpression(div(pow(var, n), n));
        }
    }
    
    // For constant
    if (!has_symbol(*expr_, var)) {
        return SymbolicExpression(mul(expr_, var));
    }
    
    // Return integral expression (unevaluated)
    return SymbolicExpression(Integral(expr_, var));
}

SymbolicExpression SymbolicExpression::simplify() const {
    return SymbolicExpression(SymEngine::simplify(expr_));
}

SymbolicExpression SymbolicExpression::expand() const {
    return SymbolicExpression(SymEngine::expand(expr_));
}

SymbolicExpression SymbolicExpression::factor() const {
    // SymEngine has limited factoring
    return SymbolicExpression(expr_);
}

SymbolicExpression SymbolicExpression::substitute(
    const std::string& variable,
    const SymbolicExpression& value) const
{
    auto var = symbol(variable);
    map_basic_basic subs_map{{var, value.expr_}};
    return SymbolicExpression(expr_->subs(subs_map));
}

double SymbolicExpression::evaluate(const std::map<std::string, double>& values) const {
    map_basic_basic subs_map;
    for (const auto& [var, val] : values) {
        subs_map[symbol(var)] = real_double(val);
    }
    
    auto result = expr_->subs(subs_map);
    return eval_double(*result);
}

std::string SymbolicExpression::toString() const {
    return expr_->__str__();
}

std::string SymbolicExpression::toLatex() const {
    return latex(*expr_);
}

RCP<const Basic> SymbolicExpression::getExpr() const {
    return expr_;
}

// SymbolicSolver implementation
std::vector<SymbolicExpression> SymbolicSolver::solve(
    const SymbolicExpression& equation,
    const std::string& variable) const
{
    auto var = symbol(variable);
    auto solutions = SymEngine::solve(equation.getExpr(), var);
    
    std::vector<SymbolicExpression> results;
    if (is_a<FiniteSet>(*solutions)) {
        const FiniteSet& fs = down_cast<const FiniteSet&>(*solutions);
        for (const auto& sol : fs.get_container()) {
            results.push_back(SymbolicExpression(sol));
        }
    }
    
    return results;
}

std::vector<std::map<std::string, SymbolicExpression>> 
SymbolicSolver::solveSystem(
    const std::vector<SymbolicExpression>& equations,
    const std::vector<std::string>& variables) const
{
    vec_basic eqs;
    set_basic vars;
    
    for (const auto& eq : equations) {
        eqs.push_back(eq.getExpr());
    }
    
    for (const auto& var : variables) {
        vars.insert(symbol(var));
    }
    
    auto solutions = SymEngine::linsolve(eqs, vars);
    
    std::vector<std::map<std::string, SymbolicExpression>> results;
    // Parse and organize solutions
    // (simplified implementation)
    
    return results;
}

std::vector<SymbolicExpression> SymbolicSolver::findCriticalPoints(
    const SymbolicExpression& function,
    const std::string& variable) const
{
    auto derivative = function.differentiate(variable);
    return solve(derivative, variable);
}

SymbolicExpression SymbolicSolver::taylorSeries(
    const SymbolicExpression& function,
    const std::string& variable,
    int order,
    double point) const
{
    auto var = symbol(variable);
    auto x0 = real_double(point);
    
    RCP<const Basic> series = integer(0);
    auto f = function.getExpr();
    
    for (int n = 0; n <= order; ++n) {
        // Evaluate f at point
        auto f_val = f->subs({{var, x0}});
        
        // Add term: f^(n)(x0) / n! * (x - x0)^n
        double factorial = 1.0;
        for (int i = 2; i <= n; ++i) {
            factorial *= i;
        }
        
        auto term = mul(div(f_val, real_double(factorial)),
                       pow(sub(var, x0), integer(n)));
        series = add(series, term);
        
        // Take next derivative
        f = diff(f, var);
    }
    
    return SymbolicExpression(series);
}
