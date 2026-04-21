#pragma once

#include <QWidget>

class QComboBox;
class QLineEdit;
class QPushButton;
class QTabWidget;
class QTextEdit;

namespace CalculatorAdvanced {

class CalculatorWidget final : public QWidget {
public:
    explicit CalculatorWidget(QWidget* parent = nullptr);
    ~CalculatorWidget() override;

private:
    void buildUi();
    void solveCurrentTab();
    void solveFunctionalTab();
    void solvePolynomialTab();
    void solveUqffTab();
    void clearAll();

    QTabWidget* tab_widget_;

    QTextEdit* functional_expression_edit_;
    QLineEdit* functional_start_edit_;
    QLineEdit* functional_end_edit_;
    QLineEdit* functional_samples_edit_;

    QTextEdit* polynomial_coefficients_edit_;

    QComboBox* uqff_equation_combo_;
    QTextEdit* uqff_parameters_edit_;

    QTextEdit* output_edit_;
    QPushButton* solve_button_;
    QPushButton* clear_button_;
};

} // namespace CalculatorAdvanced
