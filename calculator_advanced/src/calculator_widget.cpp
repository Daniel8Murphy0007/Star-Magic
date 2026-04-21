#include "../include/calculator_widget.h"

#include <QComboBox>
#include <QFont>
#include <QFormLayout>
#include <QHBoxLayout>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QRegularExpression>
#include <QTabWidget>
#include <QTextEdit>
#include <QVBoxLayout>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cctype>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kC = 299792458.0;
constexpr double kG = 6.67430e-11;
constexpr double kMu0 = 1.25663706212e-6;

std::string trimCopy(const std::string& value) {
    const auto first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) {
        return {};
    }
    const auto last = value.find_last_not_of(" \t\r\n");
    return value.substr(first, last - first + 1);
}

QString formatDouble(double value) {
    return QString::number(value, 'g', 12);
}

QString formatComplex(const std::complex<double>& value) {
    const double real = std::abs(value.real()) < 1e-12 ? 0.0 : value.real();
    const double imag = std::abs(value.imag()) < 1e-12 ? 0.0 : value.imag();

    if (imag == 0.0) {
        return formatDouble(real);
    }
    if (real == 0.0) {
        return QString("%1i").arg(formatDouble(imag));
    }
    return QString("%1 %2 %3i")
        .arg(formatDouble(real))
        .arg(imag < 0.0 ? "-" : "+")
        .arg(formatDouble(std::abs(imag)));
}

class ExpressionParser {
public:
    ExpressionParser(std::string expression, double variable)
        : expression_(std::move(expression)), variable_(variable) {}

    double parse() {
        position_ = 0;
        const double value = parseExpression();
        skipWhitespace();
        if (position_ != expression_.size()) {
            throw std::runtime_error("unexpected trailing input");
        }
        return value;
    }

private:
    double parseExpression() {
        double value = parseTerm();
        while (true) {
            skipWhitespace();
            if (match('+')) {
                value += parseTerm();
            } else if (match('-')) {
                value -= parseTerm();
            } else {
                return value;
            }
        }
    }

    double parseTerm() {
        double value = parsePower();
        while (true) {
            skipWhitespace();
            if (match('*')) {
                value *= parsePower();
            } else if (match('/')) {
                const double divisor = parsePower();
                if (std::abs(divisor) < 1e-18) {
                    throw std::runtime_error("division by zero");
                }
                value /= divisor;
            } else {
                return value;
            }
        }
    }

    double parsePower() {
        double value = parseUnary();
        skipWhitespace();
        if (match('^')) {
            value = std::pow(value, parsePower());
        }
        return value;
    }

    double parseUnary() {
        skipWhitespace();
        if (match('+')) {
            return parseUnary();
        }
        if (match('-')) {
            return -parseUnary();
        }
        return parsePrimary();
    }

    double parsePrimary() {
        skipWhitespace();
        if (match('(')) {
            const double value = parseExpression();
            if (!match(')')) {
                throw std::runtime_error("missing closing parenthesis");
            }
            return value;
        }
        if (std::isdigit(static_cast<unsigned char>(current())) || current() == '.') {
            return parseNumber();
        }
        if (std::isalpha(static_cast<unsigned char>(current()))) {
            const std::string identifier = parseIdentifier();
            skipWhitespace();
            if (match('(')) {
                const double argument = parseExpression();
                if (!match(')')) {
                    throw std::runtime_error("missing closing parenthesis after function");
                }
                return applyFunction(identifier, argument);
            }
            if (identifier == "x") {
                return variable_;
            }
            if (identifier == "pi") {
                return kPi;
            }
            if (identifier == "e") {
                return std::exp(1.0);
            }
            throw std::runtime_error("unknown identifier: " + identifier);
        }
        throw std::runtime_error("expected number, variable, or function");
    }

    double parseNumber() {
        const std::size_t start = position_;
        while (position_ < expression_.size() &&
               (std::isdigit(static_cast<unsigned char>(expression_[position_])) ||
                expression_[position_] == '.' || expression_[position_] == 'e' ||
                expression_[position_] == 'E' || expression_[position_] == '+' ||
                expression_[position_] == '-')) {
            if ((expression_[position_] == '+' || expression_[position_] == '-') && position_ != start) {
                const char previous = expression_[position_ - 1];
                if (previous != 'e' && previous != 'E') {
                    break;
                }
            }
            ++position_;
        }
        return std::stod(expression_.substr(start, position_ - start));
    }

    std::string parseIdentifier() {
        const std::size_t start = position_;
        while (position_ < expression_.size() &&
               (std::isalnum(static_cast<unsigned char>(expression_[position_])) ||
                expression_[position_] == '_')) {
            ++position_;
        }
        return expression_.substr(start, position_ - start);
    }

    double applyFunction(const std::string& name, double argument) const {
        if (name == "sin") {
            return std::sin(argument);
        }
        if (name == "cos") {
            return std::cos(argument);
        }
        if (name == "tan") {
            return std::tan(argument);
        }
        if (name == "exp") {
            return std::exp(argument);
        }
        if (name == "log") {
            return std::log(argument);
        }
        if (name == "sqrt") {
            if (argument < 0.0) {
                throw std::runtime_error("sqrt requires a non-negative argument");
            }
            return std::sqrt(argument);
        }
        if (name == "abs") {
            return std::abs(argument);
        }
        throw std::runtime_error("unknown function: " + name);
    }

    bool match(char expected) {
        skipWhitespace();
        if (current() != expected) {
            return false;
        }
        ++position_;
        return true;
    }

    void skipWhitespace() {
        while (position_ < expression_.size() &&
               std::isspace(static_cast<unsigned char>(expression_[position_]))) {
            ++position_;
        }
    }

    char current() const {
        if (position_ >= expression_.size()) {
            return '\0';
        }
        return expression_[position_];
    }

    std::string expression_;
    double variable_;
    std::size_t position_ = 0;
};

double evaluateExpression(const std::string& expression, double x) {
    ExpressionParser parser(expression, x);
    return parser.parse();
}

std::vector<double> parseCoefficientList(const QString& input) {
    QString normalized = input;
    normalized.replace(',', ' ');
    normalized.replace(';', ' ');
    normalized.replace('\n', ' ');

    const QStringList parts = normalized.split(QRegularExpression("\\s+"), Qt::SkipEmptyParts);
    std::vector<double> coefficients;
    coefficients.reserve(static_cast<std::size_t>(parts.size()));

    for (const QString& part : parts) {
        bool ok = false;
        const double value = part.toDouble(&ok);
        if (!ok) {
            throw std::runtime_error("polynomial coefficients must be numeric");
        }
        coefficients.push_back(value);
    }

    while (coefficients.size() > 1 && std::abs(coefficients.front()) < 1e-18) {
        coefficients.erase(coefficients.begin());
    }
    if (coefficients.size() < 2) {
        throw std::runtime_error("enter at least two coefficients");
    }
    return coefficients;
}

std::complex<double> evaluatePolynomial(const std::vector<double>& coefficients,
                                        const std::complex<double>& x) {
    std::complex<double> result = 0.0;
    for (double coefficient : coefficients) {
        result = result * x + coefficient;
    }
    return result;
}

std::vector<std::complex<double>> solvePolynomial(const std::vector<double>& coefficients) {
    const std::size_t degree = coefficients.size() - 1;
    std::vector<std::complex<double>> roots;
    roots.reserve(degree);

    const double radius = 1.0 + std::abs(coefficients.back() / coefficients.front());
    for (std::size_t index = 0; index < degree; ++index) {
        const double angle = 2.0 * kPi * static_cast<double>(index) / static_cast<double>(degree);
        roots.emplace_back(std::polar(radius, angle));
    }

    for (int iteration = 0; iteration < 200; ++iteration) {
        double max_delta = 0.0;
        for (std::size_t i = 0; i < degree; ++i) {
            std::complex<double> denominator = 1.0;
            for (std::size_t j = 0; j < degree; ++j) {
                if (i != j) {
                    denominator *= (roots[i] - roots[j]);
                }
            }
            if (std::abs(denominator) < 1e-18) {
                denominator = std::complex<double>(1e-9, 1e-9);
            }
            const std::complex<double> correction = evaluatePolynomial(coefficients, roots[i]) / denominator;
            roots[i] -= correction;
            max_delta = std::max(max_delta, std::abs(correction));
        }
        if (max_delta < 1e-12) {
            break;
        }
    }

    std::sort(roots.begin(), roots.end(), [](const auto& left, const auto& right) {
        if (std::abs(left.real() - right.real()) > 1e-9) {
            return left.real() < right.real();
        }
        return left.imag() < right.imag();
    });
    return roots;
}

std::map<std::string, double> parseKeyValueParameters(const QString& input) {
    QString normalized = input;
    normalized.replace('{', ' ');
    normalized.replace('}', ' ');
    normalized.replace('"', ' ');

    std::map<std::string, double> values;
    const QStringList lines = normalized.split('\n', Qt::SkipEmptyParts);
    for (QString line : lines) {
        line = line.trimmed();
        if (line.isEmpty()) {
            continue;
        }
        line.replace(',', ' ');
        line.replace(':', '=');
        const int separator = line.indexOf('=');
        if (separator <= 0) {
            continue;
        }
        const QString key = line.left(separator).trimmed();
        const QString raw_value = line.mid(separator + 1).trimmed();
        bool ok = false;
        const double value = raw_value.toDouble(&ok);
        if (ok) {
            values.emplace(key.toStdString(), value);
        }
    }
    return values;
}

double readParameter(const std::map<std::string, double>& values,
                     const std::string& key,
                     double default_value) {
    const auto found = values.find(key);
    return found == values.end() ? default_value : found->second;
}

QString solveUqffEquation(const QString& equation_name, const std::map<std::string, double>& values) {
    const double mass = readParameter(values, "M", 1.989e30);
    const double radius = std::max(readParameter(values, "r", 1.0e9), 1e-18);
    const double magnetic_field = readParameter(values, "B", 1.0e4);
    const double density = std::max(readParameter(values, "rho", 1.0), 1e-18);
    const double alpha = readParameter(values, "alpha", 0.1);
    const double hubble = readParameter(values, "H", 2.2e-18);
    const double beta = readParameter(values, "beta", 1e-12);
    const double spin = readParameter(values, "spin", 0.5);
    const double particle_mass = std::max(readParameter(values, "m", 9.109e-31), 1e-18);
    const double charge = std::max(std::abs(readParameter(values, "q", 1.602e-19)), 1e-18);
    const double lambda = std::max(readParameter(values, "lambda", 1e-9), 0.0);
    const double output_power = readParameter(values, "output_power", 1.0);
    const double input_power = std::max(readParameter(values, "input_power", 1.0), 1e-18);

    double result = 0.0;
    QString formula;

    if (equation_name == "F_U_Bi_i") {
        result = (kG * mass) / (radius * radius) + (alpha * magnetic_field) / density;
        formula = "F_U_Bi_i = G*M/r^2 + alpha*B/rho";
    } else if (equation_name == "Um") {
        result = (magnetic_field * magnetic_field) / (2.0 * kMu0);
        formula = "Um = B^2 / (2*mu0)";
    } else if (equation_name == "g_MUGE_H") {
        result = (kG * mass) / (radius * radius) * (1.0 + alpha * hubble * radius / kC);
        formula = "g_MUGE_H = (G*M/r^2) * (1 + alpha*H*r/c)";
    } else if (equation_name == "g_Magnetar") {
        result = (kG * mass) / (radius * radius) + beta * magnetic_field;
        formula = "g_Magnetar = G*M/r^2 + beta*B";
    } else if (equation_name == "g_SgrA") {
        result = (kG * mass) / (radius * radius) * (1.0 + spin * spin);
        formula = "g_SgrA = (G*M/r^2) * (1 + spin^2)";
    } else if (equation_name == "P_alpha") {
        const double omega = std::max(readParameter(values, "omega", 1.0), 1e-18);
        result = 2.0 * kPi / omega;
        formula = "P_alpha = 2*pi/omega";
    } else if (equation_name == "R_EU") {
        result = std::sqrt((2.0 * kG * mass) / (kC * kC));
        formula = "R_EU = sqrt(2*G*M/c^2)";
    } else if (equation_name == "tau_gyro") {
        result = (2.0 * kPi * particle_mass) /
                 (charge * std::max(std::abs(magnetic_field), 1e-18));
        formula = "tau_gyro = 2*pi*m/(q*B)";
    } else if (equation_name == "g_compressed") {
        result = (kG * mass) / (radius * radius) * std::exp(-lambda * radius);
        formula = "g_compressed = (G*M/r^2) * exp(-lambda*r)";
    } else if (equation_name == "eta_LENR") {
        result = output_power / input_power;
        formula = "eta_LENR = output_power / input_power";
    } else {
        throw std::runtime_error("unknown UQFF equation");
    }

    QString output;
    output += "UQFF equation\n";
    output += "-------------\n";
    output += QString("Name: %1\n").arg(equation_name);
    output += QString("Formula: %1\n").arg(formula);
    output += QString("Value: %1\n\n").arg(formatDouble(result));
    output += "Parameters used\n";
    output += "---------------\n";
    for (const auto& [key, value] : values) {
        output += QString("%1 = %2\n").arg(QString::fromStdString(key), formatDouble(value));
    }
    return output;
}

} // namespace

namespace CalculatorAdvanced {

CalculatorWidget::CalculatorWidget(QWidget* parent)
    : QWidget(parent),
      tab_widget_(nullptr),
      functional_expression_edit_(nullptr),
      functional_start_edit_(nullptr),
      functional_end_edit_(nullptr),
      functional_samples_edit_(nullptr),
      polynomial_coefficients_edit_(nullptr),
      uqff_equation_combo_(nullptr),
      uqff_parameters_edit_(nullptr),
      output_edit_(nullptr),
      solve_button_(nullptr),
      clear_button_(nullptr) {
    buildUi();
}

CalculatorWidget::~CalculatorWidget() = default;

void CalculatorWidget::buildUi() {
    auto* main_layout = new QVBoxLayout(this);

    auto* title = new QLabel("Advanced UQFF Calculator", this);
    QFont title_font = title->font();
    title_font.setPointSize(title_font.pointSize() + 2);
    title_font.setBold(true);
    title->setFont(title_font);
    main_layout->addWidget(title);

    auto* subtitle = new QLabel(
        "Functional sampling, polynomial roots, and core UQFF equations in one tab.", this);
    subtitle->setWordWrap(true);
    main_layout->addWidget(subtitle);

    tab_widget_ = new QTabWidget(this);

    auto* functional_tab = new QWidget(this);
    auto* functional_layout = new QVBoxLayout(functional_tab);
    functional_expression_edit_ = new QTextEdit(functional_tab);
    functional_expression_edit_->setPlaceholderText(
        "Enter an expression in x, for example:\n"
        "sin(x) + x^2 / 4\n"
        "exp(-x^2) * cos(x)");
    functional_expression_edit_->setMaximumHeight(110);
    functional_layout->addWidget(functional_expression_edit_);

    auto* functional_form = new QFormLayout();
    functional_start_edit_ = new QLineEdit("-5", functional_tab);
    functional_end_edit_ = new QLineEdit("5", functional_tab);
    functional_samples_edit_ = new QLineEdit("11", functional_tab);
    functional_form->addRow("Range start", functional_start_edit_);
    functional_form->addRow("Range end", functional_end_edit_);
    functional_form->addRow("Samples", functional_samples_edit_);
    functional_layout->addLayout(functional_form);
    functional_layout->addStretch();
    tab_widget_->addTab(functional_tab, "Functional");

    auto* polynomial_tab = new QWidget(this);
    auto* polynomial_layout = new QVBoxLayout(polynomial_tab);
    auto* polynomial_label = new QLabel(
        "Enter coefficients from highest degree to constant term, separated by commas or spaces.",
        polynomial_tab);
    polynomial_label->setWordWrap(true);
    polynomial_layout->addWidget(polynomial_label);
    polynomial_coefficients_edit_ = new QTextEdit(polynomial_tab);
    polynomial_coefficients_edit_->setPlaceholderText("1, -6, 11, -6");
    polynomial_coefficients_edit_->setMaximumHeight(90);
    polynomial_layout->addWidget(polynomial_coefficients_edit_);
    polynomial_layout->addStretch();
    tab_widget_->addTab(polynomial_tab, "Polynomial");

    auto* uqff_tab = new QWidget(this);
    auto* uqff_layout = new QVBoxLayout(uqff_tab);
    uqff_equation_combo_ = new QComboBox(uqff_tab);
    uqff_equation_combo_->addItems({
        "F_U_Bi_i",
        "Um",
        "g_MUGE_H",
        "g_Magnetar",
        "g_SgrA",
        "P_alpha",
        "R_EU",
        "tau_gyro",
        "g_compressed",
        "eta_LENR",
    });
    uqff_layout->addWidget(uqff_equation_combo_);

    uqff_parameters_edit_ = new QTextEdit(uqff_tab);
    uqff_parameters_edit_->setPlaceholderText(
        "M = 1.989e30\n"
        "r = 7.0e8\n"
        "B = 1.0e5\n"
        "alpha = 0.2");
    uqff_parameters_edit_->setMaximumHeight(130);
    uqff_layout->addWidget(uqff_parameters_edit_);
    uqff_layout->addStretch();
    tab_widget_->addTab(uqff_tab, "UQFF");

    main_layout->addWidget(tab_widget_);

    output_edit_ = new QTextEdit(this);
    output_edit_->setReadOnly(true);
    output_edit_->setMinimumHeight(250);
    main_layout->addWidget(output_edit_);

    auto* button_layout = new QHBoxLayout();
    solve_button_ = new QPushButton("Solve", this);
    clear_button_ = new QPushButton("Clear", this);
    button_layout->addWidget(solve_button_);
    button_layout->addWidget(clear_button_);
    button_layout->addStretch();
    main_layout->addLayout(button_layout);

    connect(solve_button_, &QPushButton::clicked, this, [this]() { solveCurrentTab(); });
    connect(clear_button_, &QPushButton::clicked, this, [this]() { clearAll(); });
}

void CalculatorWidget::solveCurrentTab() {
    try {
        switch (tab_widget_->currentIndex()) {
            case 0:
                solveFunctionalTab();
                break;
            case 1:
                solvePolynomialTab();
                break;
            case 2:
                solveUqffTab();
                break;
            default:
                output_edit_->setPlainText("Unsupported tab.");
                break;
        }
    } catch (const std::exception& error) {
        output_edit_->setPlainText(QString("Error: %1").arg(error.what()));
    }
}

void CalculatorWidget::solveFunctionalTab() {
    const std::string expression = trimCopy(functional_expression_edit_->toPlainText().toStdString());
    if (expression.empty()) {
        throw std::runtime_error("enter an expression before solving");
    }

    bool ok_start = false;
    bool ok_end = false;
    bool ok_samples = false;
    const double start = functional_start_edit_->text().toDouble(&ok_start);
    const double end = functional_end_edit_->text().toDouble(&ok_end);
    const int samples = functional_samples_edit_->text().toInt(&ok_samples);
    if (!ok_start || !ok_end || !ok_samples || samples < 2) {
        throw std::runtime_error("functional range requires valid numeric inputs and at least 2 samples");
    }

    const double step = (end - start) / static_cast<double>(samples - 1);
    const double midpoint = (start + end) * 0.5;
    const double delta = std::max(std::abs(step) * 0.5, 1e-5);

    QString output;
    output += QString("Expression: %1\n").arg(QString::fromStdString(expression));
    output += QString("Range: [%1, %2]\n").arg(formatDouble(start), formatDouble(end));
    output += QString("Samples: %1\n\n").arg(samples);
    output += "x\tf(x)\n";
    output += "-------------------------\n";

    double min_value = std::numeric_limits<double>::infinity();
    double max_value = -std::numeric_limits<double>::infinity();
    double integral = 0.0;
    double previous_value = 0.0;

    for (int index = 0; index < samples; ++index) {
        const double x = start + step * static_cast<double>(index);
        const double value = evaluateExpression(expression, x);
        output += QString("%1\t%2\n").arg(formatDouble(x), formatDouble(value));
        min_value = std::min(min_value, value);
        max_value = std::max(max_value, value);
        if (index > 0) {
            integral += 0.5 * (previous_value + value) * step;
        }
        previous_value = value;
    }

    const double derivative =
        (evaluateExpression(expression, midpoint + delta) - evaluateExpression(expression, midpoint - delta)) /
        (2.0 * delta);

    output += "\nSummary\n";
    output += "-------\n";
    output += QString("Approx. derivative at midpoint: %1\n").arg(formatDouble(derivative));
    output += QString("Approx. integral over range: %1\n").arg(formatDouble(integral));
    output += QString("Minimum sampled value: %1\n").arg(formatDouble(min_value));
    output += QString("Maximum sampled value: %1\n").arg(formatDouble(max_value));

    output_edit_->setPlainText(output);
}

void CalculatorWidget::solvePolynomialTab() {
    const std::vector<double> coefficients = parseCoefficientList(polynomial_coefficients_edit_->toPlainText());
    const std::vector<std::complex<double>> roots = solvePolynomial(coefficients);

    QString output;
    output += QString("Polynomial degree: %1\n").arg(static_cast<int>(coefficients.size() - 1));
    output += "Coefficients\n";
    output += "------------\n";
    for (std::size_t index = 0; index < coefficients.size(); ++index) {
        output += QString("a%1 = %2\n")
                      .arg(static_cast<int>(coefficients.size() - index - 1))
                      .arg(formatDouble(coefficients[index]));
    }
    output += "\nRoots\n";
    output += "-----\n";
    for (std::size_t index = 0; index < roots.size(); ++index) {
        const auto residual = evaluatePolynomial(coefficients, roots[index]);
        output += QString("r%1 = %2\n")
                      .arg(static_cast<int>(index + 1))
                      .arg(formatComplex(roots[index]));
        output += QString("  residual = %1\n").arg(formatComplex(residual));
    }

    output_edit_->setPlainText(output);
}

void CalculatorWidget::solveUqffTab() {
    const auto values = parseKeyValueParameters(uqff_parameters_edit_->toPlainText());
    output_edit_->setPlainText(solveUqffEquation(uqff_equation_combo_->currentText(), values));
}

void CalculatorWidget::clearAll() {
    functional_expression_edit_->clear();
    polynomial_coefficients_edit_->clear();
    uqff_parameters_edit_->clear();
    output_edit_->clear();
}

} // namespace CalculatorAdvanced
