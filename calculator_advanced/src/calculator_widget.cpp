#include "../include/calculator_widget.h"
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QGroupBox>
#include <QPushButton>
#include <QLineEdit>
#include <QComboBox>
#include <QWebEngineView>
#include <QLabel>
#include <QMessageBox>

CalculatorWidget::CalculatorWidget(QWidget* parent)
    : QWidget(parent)
{
    setupUI();
    initializeConnections();
}

CalculatorWidget::~CalculatorWidget() = default;

void CalculatorWidget::setupUI() {
    QVBoxLayout* mainLayout = new QVBoxLayout(this);
    
    // Title
    QLabel* title = new QLabel("Advanced UQFF Calculator", this);
    QFont titleFont;
    titleFont.setPointSize(16);
    titleFont.setBold(true);
    title->setFont(titleFont);
    title->setAlignment(Qt::AlignCenter);
    mainLayout->addWidget(title);
    
    // Create tab widget
    tabWidget_ = new QTabWidget(this);
    
    // Tab 1: Functional Equations
    QWidget* functionalTab = createFunctionalTab();
    tabWidget_->addTab(functionalTab, "Functional");
    
    // Tab 2: Polynomial Equations
    QWidget* polynomialTab = createPolynomialTab();
    tabWidget_->addTab(polynomialTab, "Polynomial");
    
    // Tab 3: UQFF Equations
    QWidget* uqffTab = createUQFFTab();
    tabWidget_->addTab(uqffTab, "UQFF");
    
    // Tab 4: Series & Derivatives
    QWidget* calculusTab = createCalculusTab();
    tabWidget_->addTab(calculusTab, "Calculus");
    
    // Tab 5: ODE & Parametric
    QWidget* odeTab = createODETab();
    tabWidget_->addTab(odeTab, "ODE/Parametric");
    
    // Tab 6: Series Expansions
    QWidget* seriesTab = createSeriesTab();
    tabWidget_->addTab(seriesTab, "Series");
    
    mainLayout->addWidget(tabWidget_);
    
    // Output area
    QGroupBox* outputGroup = new QGroupBox("Results", this);
    QVBoxLayout* outputLayout = new QVBoxLayout(outputGroup);
    
    outputView_ = new QWebEngineView(this);
    outputView_->setMinimumHeight(300);
    outputLayout->addWidget(outputView_);
    
    mainLayout->addWidget(outputGroup);
    
    // Button panel
    QHBoxLayout* buttonLayout = new QHBoxLayout();
    
    solveButton_ = new QPushButton("Solve", this);
    solveButton_->setStyleSheet("QPushButton { background-color: #4CAF50; color: white; font-weight: bold; padding: 10px; }");
    
    clearButton_ = new QPushButton("Clear", this);
    exportButton_ = new QPushButton("Export", this);
    
    buttonLayout->addWidget(solveButton_);
    buttonLayout->addWidget(clearButton_);
    buttonLayout->addWidget(exportButton_);
    buttonLayout->addStretch();
    
    mainLayout->addLayout(buttonLayout);
}

QWidget* CalculatorWidget::createFunctionalTab() {
    QWidget* tab = new QWidget();
    QVBoxLayout* layout = new QVBoxLayout(tab);
    
    // Input field
    QLabel* label = new QLabel("Enter equation (e.g., f(x) = x^2 + 2*x + 1):", tab);
    layout->addWidget(label);
    
    functionalInput_ = new QTextEdit(tab);
    functionalInput_->setMaximumHeight(100);
    functionalInput_->setPlaceholderText("f(x) = x^2 + 2*x + 1\nvar range: x = -10 to 10");
    layout->addWidget(functionalInput_);
    
    // Helper buttons
    QHBoxLayout* helpLayout = new QHBoxLayout();
    QPushButton* sinBtn = new QPushButton("sin(x)", tab);
    QPushButton* cosBtn = new QPushButton("cos(x)", tab);
    QPushButton* expBtn = new QPushButton("exp(x)", tab);
    QPushButton* logBtn = new QPushButton("log(x)", tab);
    
    connect(sinBtn, &QPushButton::clicked, [this]() { functionalInput_->insertPlainText("sin(x)"); });
    connect(cosBtn, &QPushButton::clicked, [this]() { functionalInput_->insertPlainText("cos(x)"); });
    connect(expBtn, &QPushButton::clicked, [this]() { functionalInput_->insertPlainText("exp(x)"); });
    connect(logBtn, &QPushButton::clicked, [this]() { functionalInput_->insertPlainText("log(x)"); });
    
    helpLayout->addWidget(sinBtn);
    helpLayout->addWidget(cosBtn);
    helpLayout->addWidget(expBtn);
    helpLayout->addWidget(logBtn);
    helpLayout->addStretch();
    
    layout->addLayout(helpLayout);
    layout->addStretch();
    
    return tab;
}

QWidget* CalculatorWidget::createPolynomialTab() {
    QWidget* tab = new QWidget();
    QVBoxLayout* layout = new QVBoxLayout(tab);
    
    QLabel* label = new QLabel("Enter polynomial (degree ≤ 26):", tab);
    layout->addWidget(label);
    
    polynomialInput_ = new QTextEdit(tab);
    polynomialInput_->setMaximumHeight(100);
    polynomialInput_->setPlaceholderText("x^3 - 6*x^2 + 11*x - 6 = 0");
    layout->addWidget(polynomialInput_);
    
    layout->addStretch();
    
    return tab;
}

QWidget* CalculatorWidget::createUQFFTab() {
    QWidget* tab = new QWidget();
    QVBoxLayout* layout = new QVBoxLayout(tab);
    
    // Equation selector
    QLabel* label = new QLabel("Select UQFF Equation:", tab);
    layout->addWidget(label);
    
    uqffSelector_ = new QComboBox(tab);
    
    // Populate with equations from catalog
    auto equations = uqffCatalog_.listEquations();
    for (const auto& eq : equations) {
        uqffSelector_->addItem(QString::fromStdString(eq));
    }
    
    layout->addWidget(uqffSelector_);
    
    // Parameters input
    QLabel* paramLabel = new QLabel("Parameters (JSON format):", tab);
    layout->addWidget(paramLabel);
    
    uqffParamsInput_ = new QTextEdit(tab);
    uqffParamsInput_->setMaximumHeight(150);
    uqffParamsInput_->setPlaceholderText("{\n  \"r\": 1e-15,\n  \"t\": 1e-21,\n  \"M\": 1.67e-27\n}");
    layout->addWidget(uqffParamsInput_);
    
    // Info display
    connect(uqffSelector_, QOverload<int>::of(&QComboBox::currentIndexChanged),
            [this](int index) {
        if (index >= 0) {
            QString eqName = uqffSelector_->currentText();
            const auto* eq = uqffCatalog_.getEquation(eqName.toStdString());
            if (eq) {
                QString info = QString("Description: %1\n\nLaTeX: %2\n\nRequired parameters: %3")
                    .arg(QString::fromStdString(eq->description))
                    .arg(QString::fromStdString(eq->latex))
                    .arg(QString::fromStdString(std::accumulate(
                        eq->parameters.begin(), eq->parameters.end(), std::string{},
                        [](const std::string& a, const std::string& b) {
                            return a.empty() ? b : a + ", " + b;
                        })));
                
                QMessageBox::information(this, "Equation Info", info);
            }
        }
    });
    
    layout->addStretch();
    
    return tab;
}

QWidget* CalculatorWidget::createCalculusTab() {
    QWidget* tab = new QWidget();
    QVBoxLayout* layout = new QVBoxLayout(tab);
    
    QLabel* label = new QLabel("Enter expression for calculus operations:", tab);
    layout->addWidget(label);
    
    calculusInput_ = new QTextEdit(tab);
    calculusInput_->setMaximumHeight(100);
    calculusInput_->setPlaceholderText("Examples:\n∫ x^2 dx\nd/dx sin(x)\n∂/∂x ∂/∂y (x^2*y)");
    layout->addWidget(calculusInput_);
    
    // Operation selector
    QHBoxLayout* opLayout = new QHBoxLayout();
    QLabel* opLabel = new QLabel("Operation:", tab);
    opLayout->addWidget(opLabel);
    
    QComboBox* opSelector = new QComboBox(tab);
    opSelector->addItems({"Differentiate", "Integrate", "Partial Derivative", "Gradient"});
    opLayout->addWidget(opSelector);
    opLayout->addStretch();
    
    layout->addLayout(opLayout);
    layout->addStretch();
    
    return tab;
}

QWidget* CalculatorWidget::createODETab() {
    QWidget* tab = new QWidget();
    QVBoxLayout* layout = new QVBoxLayout(tab);
    
    QLabel* label = new QLabel("Enter ODE or parametric equation:", tab);
    layout->addWidget(label);
    
    odeInput_ = new QTextEdit(tab);
    odeInput_->setMaximumHeight(100);
    odeInput_->setPlaceholderText("dy/dt = -k*y\nx(t) = cos(w*t)\ny(t) = sin(w*t)");
    layout->addWidget(odeInput_);
    
    // Initial conditions
    QLabel* icLabel = new QLabel("Initial conditions:", tab);
    layout->addWidget(icLabel);
    
    QLineEdit* icInput = new QLineEdit(tab);
    icInput->setPlaceholderText("y(0) = 1, t_max = 10");
    layout->addWidget(icInput);
    
    layout->addStretch();
    
    return tab;
}

QWidget* CalculatorWidget::createSeriesTab() {
    QWidget* tab = new QWidget();
    QVBoxLayout* layout = new QVBoxLayout(tab);
    
    QLabel* label = new QLabel("Taylor/Maclaurin series expansion:", tab);
    layout->addWidget(label);
    
    seriesInput_ = new QTextEdit(tab);
    seriesInput_->setMaximumHeight(80);
    seriesInput_->setPlaceholderText("series(exp(x), x, 5) for Taylor series of e^x to order 5");
    layout->addWidget(seriesInput_);
    
    QHBoxLayout* orderLayout = new QHBoxLayout();
    QLabel* orderLabel = new QLabel("Order:", tab);
    QSpinBox* orderSpinBox = new QSpinBox(tab);
    orderSpinBox->setRange(0, 20);
    orderSpinBox->setValue(5);
    
    orderLayout->addWidget(orderLabel);
    orderLayout->addWidget(orderSpinBox);
    orderLayout->addStretch();
    
    layout->addLayout(orderLayout);
    layout->addStretch();
    
    return tab;
}

void CalculatorWidget::initializeConnections() {
    connect(solveButton_, &QPushButton::clicked, this, &CalculatorWidget::onSolveClicked);
    connect(clearButton_, &QPushButton::clicked, this, &CalculatorWidget::onClearClicked);
    connect(exportButton_, &QPushButton::clicked, this, &CalculatorWidget::onExportClicked);
}

void CalculatorWidget::onSolveClicked() {
    int currentTab = tabWidget_->currentIndex();
    
    switch (currentTab) {
        case 0: // Functional
            solveFunctional();
            break;
        case 1: // Polynomial
            solvePolynomial();
            break;
        case 2: // UQFF
            solveUQFF();
            break;
        case 3: // Calculus
            solveCalculus();
            break;
        case 4: // ODE
            solveODE();
            break;
        case 5: // Series
            solveSeries();
            break;
        default:
            break;
    }
}

void CalculatorWidget::solveFunctional() {
    std::string input = functionalInput_->toPlainText().toStdString();
    
    try {
        // Parse with ANTLR4
        auto parsed = parser_.parse(input);
        
        if (!parsed.isValid) {
            displayError(parsed.errorMessage);
            return;
        }
        
        // Use SymEngine to solve
        SymbolicExpression expr(parsed.rawInput);
        auto simplified = expr.simplify();
        
        // Display result
        QString html = QString("<html><head><script type=\"text/javascript\" async "
                             "src=\"https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.5/MathJax.js?config=TeX-MML-AM_CHTML\">"
                             "</script></head><body>"
                             "<h3>Solution:</h3>"
                             "<p>\\[%1\\]</p>"
                             "</body></html>").arg(QString::fromStdString(simplified.toLatex()));
        
        outputView_->setHtml(html);
        
    } catch (const std::exception& e) {
        displayError(std::string("Error: ") + e.what());
    }
}

void CalculatorWidget::solvePolynomial() {
    std::string input = polynomialInput_->toPlainText().toStdString();
    
    try {
        // Parse and solve polynomial
        auto parsed = parser_.parse(input);
        
        if (!parsed.isValid) {
            displayError(parsed.errorMessage);
            return;
        }
        
        SymbolicExpression expr(parsed.rawInput);
        SymbolicSolver solver;
        
        // Find variable (assume 'x' for now)
        auto solutions = solver.solve(expr, "x");
        
        // Display roots
        QString html = "<html><body><h3>Polynomial Roots:</h3><ul>";
        for (const auto& sol : solutions) {
            html += QString("<li>\\(%1\\)</li>").arg(QString::fromStdString(sol.toLatex()));
        }
        html += "</ul></body></html>";
        
        outputView_->setHtml(html);
        
    } catch (const std::exception& e) {
        displayError(std::string("Error: ") + e.what());
    }
}

void CalculatorWidget::solveUQFF() {
    try {
        std::string eqName = uqffSelector_->currentText().toStdString();
        std::string paramsText = uqffParamsInput_->toPlainText().toStdString();
        
        // Parse JSON parameters (simplified - should use proper JSON parser)
        std::map<std::string, double> params;
        // ... parse JSON ...
        
        double result = uqffCatalog_.calculate(eqName, params);
        
        QString html = QString("<html><body><h3>UQFF Calculation Result:</h3>"
                             "<p><strong>Equation:</strong> %1</p>"
                             "<p><strong>Result:</strong> %2</p>"
                             "</body></html>")
            .arg(QString::fromStdString(eqName))
            .arg(result, 0, 'e', 6);
        
        outputView_->setHtml(html);
        
    } catch (const std::exception& e) {
        displayError(std::string("Error: ") + e.what());
    }
}

void CalculatorWidget::solveCalculus() {
    std::string input = calculusInput_->toPlainText().toStdString();
    
    try {
        auto parsed = parser_.parse(input);
        
        if (!parsed.isValid) {
            displayError(parsed.errorMessage);
            return;
        }
        
        SymbolicExpression expr(parsed.rawInput);
        
        // Determine operation based on input
        SymbolicExpression result = expr;
        if (parsed.type == EquationType::DERIVATIVE) {
            result = expr.differentiate("x"); // Assume x for simplicity
        } else if (parsed.type == EquationType::INTEGRAL) {
            result = expr.integrate("x");
        }
        
        QString html = QString("<html><body><h3>Result:</h3>"
                             "<p>\\[%1\\]</p>"
                             "</body></html>").arg(QString::fromStdString(result.toLatex()));
        
        outputView_->setHtml(html);
        
    } catch (const std::exception& e) {
        displayError(std::string("Error: ") + e.what());
    }
}

void CalculatorWidget::solveODE() {
    displayError("ODE solving not yet implemented - requires numerical integration");
}

void CalculatorWidget::solveSeries() {
    displayError("Series expansion not yet implemented");
}

void CalculatorWidget::onClearClicked() {
    functionalInput_->clear();
    polynomialInput_->clear();
    uqffParamsInput_->clear();
    calculusInput_->clear();
    odeInput_->clear();
    seriesInput_->clear();
    outputView_->setHtml("");
}

void CalculatorWidget::onExportClicked() {
    // TODO: Implement export functionality
    QMessageBox::information(this, "Export", "Export functionality coming soon!");
}

void CalculatorWidget::displayError(const std::string& message) {
    QString html = QString("<html><body><h3 style=\"color: red;\">Error</h3>"
                         "<p>%1</p>"
                         "</body></html>").arg(QString::fromStdString(message));
    outputView_->setHtml(html);
}
