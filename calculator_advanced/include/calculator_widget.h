/**
 * Calculator Widget (Main GUI)
 * Thread b6d9bc22 Priority 2 - Iteration #30-32 Complete Implementation
 * 
 * Main interactive calculator widget for source2.cpp Tab 22
 * Integrates ANTLR4, SymEngine, GSL, UQFF equations, and plotting
 */

#pragma once

#include "antlr4_parser.h"
#include "symengine_wrapper.h"
#include "polynomial_solver.h"
#include "uqff_equations.h"
#include "dimensional_analysis.h"
#include "equation_solver.h"
#include "plotter_widget.h"

#include <QWidget>
#include <QTabWidget>
#include <QLineEdit>
#include <QPushButton>
#include <QTextEdit>
#include <QComboBox>
#include <QWebEngineView>
#include <memory>

namespace CalculatorAdvanced {

/**
 * @brief Main calculator GUI widget
 */
class CalculatorWidget : public QWidget {
    Q_OBJECT
    
public:
    explicit CalculatorWidget(QWidget* parent = nullptr);
    ~CalculatorWidget() override;
    
    /**
     * @brief Initialize calculator with default settings
     */
    void initialize();
    
public slots:
    /**
     * @brief Solve button clicked
     */
    void onSolveClicked();
    
    /**
     * @brief Clear input/output
     */
    void onClearClicked();
    
    /**
     * @brief UQFF equation selected from catalog
     */
    void onUQFFEquationSelected(const QString& equation_name);
    
    /**
     * @brief Plot button clicked
     */
    void onPlotClicked();
    
    /**
     * @brief Export result to file
     */
    void onExportClicked();
    
    /**
     * @brief Tab changed (switch between equation types)
     */
    void onTabChanged(int index);
    
signals:
    void solutionComputed(const QString& result);
    void errorOccurred(const QString& error);
    
private:
    // UI Components
    QTabWidget* tab_widget_;
    
    // Tab 1: Functional Equations
    QWidget* functional_tab_;
    QLineEdit* functional_input_;
    QPushButton* functional_solve_btn_;
    QTextEdit* functional_output_;
    QWebEngineView* functional_latex_view_;  // MathJax rendering
    
    // Tab 2: Parametric Curves
    QWidget* parametric_tab_;
    QLineEdit* parametric_x_input_;
    QLineEdit* parametric_y_input_;
    QLineEdit* parametric_z_input_;
    QLineEdit* parametric_t_min_;
    QLineEdit* parametric_t_max_;
    QPushButton* parametric_solve_btn_;
    PlotterWidget* parametric_plotter_;
    
    // Tab 3: ODE Systems
    QWidget* ode_tab_;
    QTextEdit* ode_input_;  // Multi-line for system of ODEs
    QTextEdit* ode_initial_conditions_;
    QLineEdit* ode_t_min_;
    QLineEdit* ode_t_max_;
    QPushButton* ode_solve_btn_;
    PlotterWidget* ode_plotter_;
    
    // Tab 4: Series Expansions
    QWidget* series_tab_;
    QLineEdit* series_function_input_;
    QLineEdit* series_variable_;
    QLineEdit* series_expansion_point_;
    QLineEdit* series_order_;
    QPushButton* series_solve_btn_;
    QTextEdit* series_output_;
    QWebEngineView* series_latex_view_;
    
    // Tab 5: Polynomial Solver
    QWidget* polynomial_tab_;
    QLineEdit* polynomial_degree_;
    QTextEdit* polynomial_coefficients_;  // Comma-separated
    QPushButton* polynomial_solve_btn_;
    QTextEdit* polynomial_output_;
    PlotterWidget* polynomial_plotter_;  // Roots on complex plane
    
    // Tab 6: UQFF Equations
    QWidget* uqff_tab_;
    QComboBox* uqff_equation_selector_;
    QTextEdit* uqff_parameters_;  // JSON-like format
    QPushButton* uqff_solve_btn_;
    QTextEdit* uqff_output_;
    QWebEngineView* uqff_latex_view_;
    PlotterWidget* uqff_plotter_;
    
    // Backend Solvers
    std::unique_ptr<EquationSolver> solver_;
    std::unique_ptr<SymbolicSolver> symbolic_solver_;
    std::unique_ptr<PolynomialSolver> polynomial_solver_;
    std::unique_ptr<UQFFEquationCatalog> uqff_catalog_;
    std::unique_ptr<DimensionalSystem> dimensional_system_;
    
    // Helper Methods
    void setupUI();
    void setupFunctionalTab();
    void setupParametricTab();
    void setupODETab();
    void setupSeriesTab();
    void setupPolynomialTab();
    void setupUQFFTab();
    
    void connectSignals();
    
    void solveFunctional();
    void solveParametric();
    void solveODE();
    void solveSeries();
    void solvePolynomial();
    void solveUQFF();
    
    void displayLatex(QWebEngineView* view, const QString& latex);
    void displaySolution(const FunctionalSolution& solution);
    void displayParametricSolution(const ParametricSolution& solution);
    void displayODESolution(const ODESolution& solution);
    void displaySeriesSolution(const SeriesSolution& solution);
    
    QString generateMathJaxHTML(const QString& latex);
};

/**
 * @brief Integration into source2.cpp Principal GUI
 * 
 * Add to source2.cpp TabWidget creation section:
 * 
 * ```cpp
 * #include "calculator_advanced/include/calculator_widget.h"
 * 
 * // In MainWindow constructor after existing tabs:
 * auto* advanced_calc = new CalculatorAdvanced::CalculatorWidget(this);
 * advanced_calc->initialize();
 * tabWidget->addTab(advanced_calc, QIcon(":/icons/calculator.png"), "⚡ Advanced Calculator");
 * ```
 */

/**
 * @brief Example usage (standalone):
 * 
 * int main(int argc, char* argv[]) {
 *     QApplication app(argc, argv);
 *     
 *     CalculatorWidget* calculator = new CalculatorWidget();
 *     calculator->initialize();
 *     calculator->resize(1200, 800);
 *     calculator->show();
 *     
 *     return app.exec();
 * }
 */

/**
 * @brief Advanced Features from Iteration #32:
 * 
 * 1. Voice Recognition (Qt Speech module):
 *    - "Solve x squared plus y squared equals 25"
 *    - "Plot F U Bi i from 1e14 to 1e16 meters"
 * 
 * 2. VR/AR Export (OpenXR integration):
 *    - Export 3D parametric curves to .gltf
 *    - Interactive equation manipulation in virtual space
 * 
 * 3. GPU Acceleration (CUDA/OpenCL):
 *    - Offload ODE integration to GPU
 *    - Parallel polynomial root finding
 * 
 * 4. Blockchain Collaboration (IPFS + Ethereum):
 *    - Share UQFF solutions with cryptographic verification
 *    - Version control for equation catalogs
 * 
 * 5. Gamification:
 *    - Achievement system for solving complex UQFF equations
 *    - Leaderboards for fastest convergence times
 * 
 * 6. Multilingual Support:
 *    - Equation parsing in LaTeX, MathML, Wolfram Language
 *    - UI translations (English, Spanish, Chinese, Japanese)
 * 
 * 7. Jupyter Integration:
 *    - Export to .ipynb with interactive widgets
 *    - Embed calculator in Jupyter notebooks
 * 
 * These features are scaffolded in Iteration #32 code blocks but require
 * additional dependencies (qt6-speech, openxr, cuda, ipfs-api).
 */

} // namespace CalculatorAdvanced
