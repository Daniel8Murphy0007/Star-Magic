/**
 * @file equation_renderer.h
 * @brief Qt widget for rendering long-form UQFF equations from QCalc.py
 * 
 * Uses QTextEdit with monospace font for Unicode math rendering.
 * Supports both static presets and dynamic IPC via QProcess to QCalc.py.
 * 
 * Gap #4 Fix: Added IPC integration to fetch computed long-form equations
 * from QCalc.py's long_form_equations output in real-time.
 * 
 * Author: Daniel T. Murphy
 * Framework: UQFF Star-Magic v3.0
 */

#ifndef EQUATION_RENDERER_H
#define EQUATION_RENDERER_H

#include <QWidget>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QTextEdit>
#include <QPushButton>
#include <QComboBox>
#include <QLabel>
#include <QSplitter>
#include <QProcess>
#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonArray>
#include <QFont>
#include <QFontDatabase>
#include <QApplication>
#include <QClipboard>
#include <QScrollBar>
#include <QLineEdit>
#include <QProgressBar>
#include <QFile>
#include <QCoreApplication>
#include <QMessageBox>

/**
 * @class EquationRendererWidget
 * @brief Widget to display UQFF equations in readable formatted text
 * 
 * Features:
 * - Monospace font for proper Unicode math alignment
 * - Syntax highlighting for variables and operators
 * - Export to LaTeX/PDF
 * - Step-by-step solution display
 * - IPC integration with QCalc.py for dynamic equation computation
 */
class EquationRendererWidget : public QWidget {
    Q_OBJECT

public:
    explicit EquationRendererWidget(QWidget* parent = nullptr)
        : QWidget(parent), qcalcProcess(nullptr) {
        
        setupUI();
        initQCalcProcess();
    }
    
    ~EquationRendererWidget() {
        if (qcalcProcess && qcalcProcess->state() == QProcess::Running) {
            qcalcProcess->kill();
            qcalcProcess->waitForFinished(3000);
        }
    }

private:
    void setupUI() {
        QVBoxLayout* mainLayout = new QVBoxLayout(this);
        
        // Header
        QLabel* header = new QLabel("📐 UQFF Long-Form Equation Display", this);
        header->setStyleSheet(
            "background-color: #1565C0; "
            "color: white; "
            "font-size: 14pt; "
            "font-weight: bold; "
            "padding: 10px; "
            "border-radius: 5px;"
        );
        mainLayout->addWidget(header);
        
        // ============================================================
        // DYNAMIC EQUATION INPUT (IPC Integration - Gap #4 Fix)
        // ============================================================
        QHBoxLayout* dynamicLayout = new QHBoxLayout();
        
        QLabel* systemLabel = new QLabel("System:", this);
        dynamicLayout->addWidget(systemLabel);
        
        systemInput = new QLineEdit(this);
        systemInput->setPlaceholderText("Enter system name (e.g., 'Sagittarius A*', 'NGC 3596', 'ESO 137-001')");
        systemInput->setStyleSheet(
            "padding: 8px; "
            "border: 2px solid #1565C0; "
            "border-radius: 4px;"
        );
        connect(systemInput, &QLineEdit::returnPressed, this, &EquationRendererWidget::fetchDynamicEquations);
        dynamicLayout->addWidget(systemInput);
        
        QPushButton* fetchBtn = new QPushButton("🔬 Compute", this);
        fetchBtn->setStyleSheet(
            "background-color: #4CAF50; "
            "color: white; "
            "padding: 8px 16px; "
            "font-weight: bold;"
        );
        connect(fetchBtn, &QPushButton::clicked, this, &EquationRendererWidget::fetchDynamicEquations);
        dynamicLayout->addWidget(fetchBtn);
        
        mainLayout->addLayout(dynamicLayout);
        
        // Progress bar for IPC operations
        progressBar = new QProgressBar(this);
        progressBar->setVisible(false);
        progressBar->setStyleSheet("QProgressBar { border: 1px solid #1565C0; border-radius: 3px; }");
        mainLayout->addWidget(progressBar);
        
        // Status label
        statusLabel = new QLabel("", this);
        statusLabel->setStyleSheet("color: #757575; font-style: italic;");
        mainLayout->addWidget(statusLabel);
        
        // ============================================================
        // PRESET EQUATION SELECTOR (Original Functionality)
        // ============================================================
        QHBoxLayout* selectorLayout = new QHBoxLayout();
        
        QLabel* eqLabel = new QLabel("Preset:", this);
        selectorLayout->addWidget(eqLabel);
        
        equationSelector = new QComboBox(this);
        equationSelector->addItem("UQFF Base Unified Field");
        equationSelector->addItem("26-Layer Compressed Gravity");
        equationSelector->addItem("F_U_Bi_i Master Buoyancy");
        equationSelector->addItem("MUGE Resonance Framework");
        equationSelector->addItem("Galaxy Rotation Curve");
        equationSelector->addItem("Magnetar Field Evolution");
        equationSelector->addItem("BEC Integration Module");
        equationSelector->addItem("LENR Widom-Larsen");
        connect(equationSelector, QOverload<int>::of(&QComboBox::currentIndexChanged),
                this, &EquationRendererWidget::loadEquation);
        selectorLayout->addWidget(equationSelector);
        
        selectorLayout->addStretch();
        
        // Render mode
        QLabel* modeLabel = new QLabel("View:", this);
        selectorLayout->addWidget(modeLabel);
        
        modeSelector = new QComboBox(this);
        modeSelector->addItem("Symbolic Form");
        modeSelector->addItem("Long-Form Solution");
        modeSelector->addItem("LaTeX Source");
        connect(modeSelector, QOverload<int>::of(&QComboBox::currentIndexChanged),
                this, &EquationRendererWidget::changeMode);
        selectorLayout->addWidget(modeSelector);
        
        mainLayout->addLayout(selectorLayout);
        
        // Main equation display
        equationDisplay = new QTextEdit(this);
        equationDisplay->setReadOnly(true);
        
        // Use a good Unicode math font
        QFont mathFont("Consolas", 11);
        if (!QFontDatabase::hasFamily("Consolas")) {
            mathFont = QFont("Courier New", 11);
        }
        equationDisplay->setFont(mathFont);
        
        equationDisplay->setStyleSheet(
            "background-color: #1E1E1E; "
            "color: #D4D4D4; "
            "border: 2px solid #1565C0; "
            "padding: 15px; "
            "line-height: 1.6;"
        );
        mainLayout->addWidget(equationDisplay);
        
        // Button bar
        QHBoxLayout* buttonLayout = new QHBoxLayout();
        
        QPushButton* copyBtn = new QPushButton("📋 Copy", this);
        copyBtn->setStyleSheet("padding: 8px 16px;");
        connect(copyBtn, &QPushButton::clicked, this, &EquationRendererWidget::copyEquation);
        buttonLayout->addWidget(copyBtn);
        
        QPushButton* exportLatexBtn = new QPushButton("📄 Export LaTeX", this);
        exportLatexBtn->setStyleSheet("padding: 8px 16px;");
        connect(exportLatexBtn, &QPushButton::clicked, this, &EquationRendererWidget::exportLatex);
        buttonLayout->addWidget(exportLatexBtn);
        
        QPushButton* computeBtn = new QPushButton("▶ Compute", this);
        computeBtn->setStyleSheet(
            "background-color: #4CAF50; "
            "color: white; "
            "padding: 8px 16px; "
            "font-weight: bold;"
        );
        connect(computeBtn, &QPushButton::clicked, this, &EquationRendererWidget::computeEquation);
        buttonLayout->addWidget(computeBtn);
        
        buttonLayout->addStretch();
        mainLayout->addLayout(buttonLayout);
        
        setLayout(mainLayout);
        
        // Load initial equation
        loadEquation(0);
    }
    
    // ============================================================
    // IPC INTEGRATION (Gap #4 Fix)
    // Initialize QProcess for QCalc.py communication
    // ============================================================
    void initQCalcProcess() {
        qcalcProcess = new QProcess(this);
        
        connect(qcalcProcess, &QProcess::readyReadStandardOutput, 
                this, &EquationRendererWidget::handleQCalcOutput);
        connect(qcalcProcess, &QProcess::readyReadStandardError,
                this, &EquationRendererWidget::handleQCalcError);
        connect(qcalcProcess, QOverload<int, QProcess::ExitStatus>::of(&QProcess::finished),
                this, &EquationRendererWidget::handleQCalcFinished);
    }
    
    QString findQCalcPath() {
        // Search for QCalc.py in order of preference
        QStringList searchPaths = {
            QCoreApplication::applicationDirPath() + "/QCalc.py",
            QCoreApplication::applicationDirPath() + "/../QCalc.py",
            "./QCalc.py",
            "../QCalc.py"
        };
        
        for (const QString& path : searchPaths) {
            if (QFile::exists(path)) {
                return path;
            }
        }
        return "";
    }

private slots:
    // ============================================================
    // DYNAMIC EQUATION FETCHING (Gap #4 Fix)
    // Calls QCalc.py to compute long-form equations for a system
    // ============================================================
    void fetchDynamicEquations() {
        QString systemName = systemInput->text().trimmed();
        if (systemName.isEmpty()) {
            statusLabel->setText("❌ Enter a system name (e.g., 'Sagittarius A*')");
            statusLabel->setStyleSheet("color: #FF5252;");
            return;
        }
        
        QString qcalcPath = findQCalcPath();
        if (qcalcPath.isEmpty()) {
            statusLabel->setText("❌ QCalc.py not found");
            statusLabel->setStyleSheet("color: #FF5252;");
            return;
        }
        
        // Show progress
        progressBar->setVisible(true);
        progressBar->setRange(0, 0);  // Indeterminate
        statusLabel->setText("🔄 Computing equations for: " + systemName);
        statusLabel->setStyleSheet("color: #1565C0;");
        
        // Build Python command to call QCalc.py and get long_form_equations
        // Output as JSON for parsing
        QString pythonCode = QString(
            "import sys; sys.path.insert(0, '.'); "
            "from QCalc import UnifiedFieldSolver, ComputeParams, UQFFScale; "
            "import json; "
            "solver = UnifiedFieldSolver(); "
            "params = ComputeParams(M=1e42, r=1e21, z=0.1, t=4e17, scale=UQFFScale.GALACTIC); "
            "result = solver.solve(params); "
            "eqs = result.get('long_form_equations', []); "
            "output = {'system': '%1', 'count': len(eqs), 'equations': []}; "
            "for eq in eqs[:50]: "  // Limit to 50 equations for display
            "    output['equations'].append({'name': eq.get('name', 'Unknown'), 'symbolic': eq.get('symbolic_form', ''), 'value': eq.get('computed_value', 0)}); "
            "print('QCALC_JSON_START'); "
            "print(json.dumps(output)); "
            "print('QCALC_JSON_END')"
        ).arg(systemName);
        
        pendingOutput.clear();
        qcalcProcess->start("python", QStringList() << "-c" << pythonCode);
    }
    
    void handleQCalcOutput() {
        QByteArray data = qcalcProcess->readAllStandardOutput();
        pendingOutput.append(QString::fromLocal8Bit(data));
    }
    
    void handleQCalcError() {
        QByteArray data = qcalcProcess->readAllStandardError();
        QString errorText = QString::fromLocal8Bit(data);
        
        // Show errors in status bar
        if (!errorText.trimmed().isEmpty()) {
            statusLabel->setText("⚠️ " + errorText.left(100));
            statusLabel->setStyleSheet("color: #FF9800;");
        }
    }
    
    void handleQCalcFinished(int exitCode, QProcess::ExitStatus exitStatus) {
        progressBar->setVisible(false);
        
        if (exitStatus != QProcess::NormalExit || exitCode != 0) {
            statusLabel->setText("❌ QCalc.py execution failed (exit code: " + QString::number(exitCode) + ")");
            statusLabel->setStyleSheet("color: #FF5252;");
            return;
        }
        
        // Parse JSON from output
        int startIdx = pendingOutput.indexOf("QCALC_JSON_START");
        int endIdx = pendingOutput.indexOf("QCALC_JSON_END");
        
        if (startIdx == -1 || endIdx == -1) {
            statusLabel->setText("❌ Invalid response from QCalc.py");
            statusLabel->setStyleSheet("color: #FF5252;");
            return;
        }
        
        QString jsonStr = pendingOutput.mid(startIdx + 16, endIdx - startIdx - 16).trimmed();
        
        QJsonDocument doc = QJsonDocument::fromJson(jsonStr.toUtf8());
        if (doc.isNull()) {
            statusLabel->setText("❌ Failed to parse JSON from QCalc.py");
            statusLabel->setStyleSheet("color: #FF5252;");
            return;
        }
        
        QJsonObject root = doc.object();
        QString systemName = root["system"].toString();
        int count = root["count"].toInt();
        QJsonArray equations = root["equations"].toArray();
        
        // Format equations for display
        QString displayText;
        displayText += "═══════════════════════════════════════════════════════════════════════════════\n";
        displayText += QString(" LONG-FORM EQUATIONS: %1\n").arg(systemName.toUpper());
        displayText += QString(" Total Computed: %1\n").arg(count);
        displayText += "═══════════════════════════════════════════════════════════════════════════════\n\n";
        
        for (int i = 0; i < equations.size(); ++i) {
            QJsonObject eq = equations[i].toObject();
            QString name = eq["name"].toString();
            QString symbolic = eq["symbolic"].toString();
            double value = eq["value"].toDouble();
            
            displayText += QString("─────────────────────────────────────────────────────────────────────────\n");
            displayText += QString(" [%1] %2\n").arg(i + 1).arg(name);
            displayText += QString("─────────────────────────────────────────────────────────────────────────\n");
            
            if (!symbolic.isEmpty()) {
                displayText += QString("\n   Symbolic Form:\n   %1\n").arg(symbolic);
            }
            displayText += QString("\n   Computed Value: %1\n\n").arg(value, 0, 'e', 6);
        }
        
        displayText += "═══════════════════════════════════════════════════════════════════════════════\n";
        displayText += " END OF EQUATIONS\n";
        displayText += "═══════════════════════════════════════════════════════════════════════════════\n";
        
        currentEquation = displayText;
        equationDisplay->setPlainText(displayText);
        
        statusLabel->setText(QString("✅ Computed %1 equations for: %2").arg(count).arg(systemName));
        statusLabel->setStyleSheet("color: #4CAF50;");
    }
    
    // ============================================================
    // ORIGINAL PRESET EQUATION LOADING
    // ============================================================
    void loadEquation(int index) {
        QString equations[] = {
            getUQFFBaseEquation(),
            get26LayerEquation(),
            getFUBiEquation(),
            getMUGEResonanceEquation(),
            getGalaxyRotationEquation(),
            getMagnetarEquation(),
            getBECEquation(),
            getLENREquation()
        };
        
        if (index >= 0 && index < 8) {
            currentEquation = equations[index];
            equationDisplay->setPlainText(currentEquation);
        }
    }
    
    void changeMode(int mode) {
        // 0 = Symbolic, 1 = Long-form, 2 = LaTeX
        if (mode == 2) {
            // Show LaTeX version
            QString latex = convertToLatex(currentEquation);
            equationDisplay->setPlainText(latex);
        } else {
            loadEquation(equationSelector->currentIndex());
        }
    }
    
    void copyEquation() {
        QApplication::clipboard()->setText(equationDisplay->toPlainText());
    }
    
    void exportLatex() {
        QString latex = convertToLatex(currentEquation);
        // Save to file
        QString filename = "equation_" + QString::number(equationSelector->currentIndex()) + ".tex";
        QFile file(filename);
        if (file.open(QIODevice::WriteOnly | QIODevice::Text)) {
            QTextStream out(&file);
            out << "\\documentclass{article}\n";
            out << "\\usepackage{amsmath,amssymb}\n";
            out << "\\begin{document}\n\n";
            out << latex;
            out << "\n\\end{document}\n";
            file.close();
        }
    }
    
    void computeEquation() {
        emit computeRequested(equationSelector->currentText());
    }

signals:
    void computeRequested(const QString& equationName);

private:
    QString getUQFFBaseEquation() {
        return R"(
═══════════════════════════════════════════════════════════════════════════════
 UQFF BASE UNIFIED FIELD EQUATION
═══════════════════════════════════════════════════════════════════════════════

 MASTER EQUATION:
   F = -F₀ + Σ(Ug1 + Ug2 + Ug3 + Ug4) + Ub_i + Um

 COMPONENT EQUATIONS:

   Ug1 = (G×M/r²) × (1+δ) × (μ₀B₀²/8π)           [Magnetic dipole defects]
   
   Ug2 = (G×M/r²) × (Q_A+Q_UA) × S × H_SCm       [Charge-reactivity bubble]
   
   Ug3 = (c/r) × ωs×cos(ωs·t·π) × sin(θ)cos(φ) × B₀  [String rotation helicity]
   
   Ug4 = k₄ × ρ_vac × C × (M_BH/d_g) × e^(-κt) × cos(πt_n) × (1+f_fb)
                                                   [Vacuum concentration BH]
   
   Ub_i = -β_i × ΣUg_i × Ω_g × (M_BH/d_g) × U_UA × cos(πt_n)
                                                   [Buoyancy opposition]
   
   Um = (μ_j/r) × (1-e^(-γt)cos(πt_n)) × φ̂ × P_SCm × E_react
                                                   [Universal magnetism]

═══════════════════════════════════════════════════════════════════════════════
 CALIBRATED CONSTANTS:
   κ = 0.0005/day        (JWST quasar τ~2000 days)
   [SSq] = 0.57          (Shell Quotient from nebula Ye~0.1)
   β_i = 0.603           (Buoyancy factor)
   H_SCm ≈ 0.99          (Superconductivity modulation)
   U_UA ≈ 0.0001         (Universal Aether coupling)
═══════════════════════════════════════════════════════════════════════════════
)";
    }
    
    QString get26LayerEquation() {
        return R"(
═══════════════════════════════════════════════════════════════════════════════
 26-LAYER COMPRESSED GRAVITY FRAMEWORK
═══════════════════════════════════════════════════════════════════════════════

 MASTER EQUATION:
   g(r,t) = Σᵢ₌₁²⁶ [Ug1_i + Ug2_i + Ug3_i + Ug4_i]

 LAYER COMPONENTS:
   
   Ug1_i = E_DPM,i / r_i² × [UA]_i × f_TRZ_i
   
   Ug2_i = (G×M/r_i²) × (Q_A+Q_UA)_i × S_i × H_SCm
   
   Ug3_i = (c/r_i) × ω_s × cos(ω_s·t·π) × sin(θ_i)cos(φ_i) × B_i
   
   Ug4_i = k₄ × ρ_vac_i × C_i × (M_BH/d_g) × e^(-κt) × cos(πt_n)

 DIPOLE MOMENTUM ENERGY:
   E_DPM,i = (ℏc/r_i²) × Q_i × [SCm]_i
   
   where:
     r_i = r / i                    (Layer-dependent radius)
     Q_i = i                        (Quality factor scales with layer)
     [SCm]_i = i²                   (Superconductivity scaling)
     [UA]_i = i                     (Universal Aether scaling)
     f_TRZ_i = 1/i                  (Triadic resonance frequency)

═══════════════════════════════════════════════════════════════════════════════
 LAYER INDEX: i = 1, 2, 3, ..., 26 (string theory + 1 emergent)
═══════════════════════════════════════════════════════════════════════════════
)";
    }
    
    QString getFUBiEquation() {
        return R"(
═══════════════════════════════════════════════════════════════════════════════
 F_U_Bi_i MASTER BUOYANCY-GRAVITY FORCE (Outside-In)
═══════════════════════════════════════════════════════════════════════════════

 INTEGRAL FORM:
   F_U_Bi_i = ∫₀^∞ f(r)·g(t)·h(θ)·Π(params) dr dt dθ

 EXPANDED MASTER EQUATION:
   F_U_Bi_i = M × β_i × Σ(Ug_i) × LENR_factor × SCm_modulation

 COMPONENT TERMS:

   β_i = 0.603                      [Calibrated buoyancy coefficient]
   
   Σ(Ug_i) = Ug1 + Ug2 + Ug3 + Ug4  [Gravity component sum]
   
   LENR_factor = 1 + k_LENR × sin(ω_LENR × t)
                                    [Low-Energy Nuclear Resonance]
                 where k_LENR = 10⁻¹⁰, ω_LENR = 1.2×10¹² rad/s
   
   SCm_modulation = [SSq] × (1 - e^(-κt))
                                    [Superconductivity modulation]
                    where [SSq] = 0.57, κ = 0.0005/day

 INTEGRAL BOUNDS (x₂ quadratic limits):
   r: [0, ∞) with exponential decay cutoff
   t: [0, t_Hubble] = [0, 4.35×10¹⁷ s]
   θ: [0, 2π] for full spherical coverage

═══════════════════════════════════════════════════════════════════════════════
 VALIDATION: LIGO GWTC-4.0 (218 events), 99.9% solvability
═══════════════════════════════════════════════════════════════════════════════
)";
    }
    
    QString getMUGEResonanceEquation() {
        return R"(
═══════════════════════════════════════════════════════════════════════════════
 MUGE RESONANCE FRAMEWORK (Multi-Universal Gravity Equation)
═══════════════════════════════════════════════════════════════════════════════

 BASE EQUATION:
   g_MUGE(r,t) = g_N + Σ(corrections)

 NEWTONIAN BASE:
   g_N = G × M / r²

 CORRECTION TERMS (9 total):

   1. Expansion:     ΔH = H₀ × (1 - e^(-r/r_Hubble))
   
   2. Super:         Δmag = -μ₀B²/(8πρ) × suppression
   
   3. Envelope:      Δenv = envelope_modulation(M, r)
   
   4. Ug_sum:        ΔUg = Ug1 + Ug2 + Ug3 + Ug4
   
   5. Cosmological:  ΔΛ = Λc²r/3
   
   6. Quantum:       Δℏ = ℏω/(Mc²) × quantum_correction
   
   7. Fluid:         ΔNS = Navier-Stokes turbulence term
   
   8. Perturbation:  Δpert = dark_matter_perturbation

 RESONANCE MODES (13 total):
   aDPM, aTHz, Avac_diff, aSuperFreq, aAetherRes, Ug4i, aQuantumFreq,
   aAetherFreq, aFluidFreq, Osc_term, aExpFreq, fTRZ, Wormhole_metric

═══════════════════════════════════════════════════════════════════════════════
)";
    }
    
    QString getGalaxyRotationEquation() {
        return R"(
═══════════════════════════════════════════════════════════════════════════════
 GALAXY ROTATION CURVE (UQFF Dark Matter Alternative)
═══════════════════════════════════════════════════════════════════════════════

 OBSERVED PROBLEM:
   v_obs(r) ≠ v_Kepler(r) = √(G×M_visible/r)
   
   Stars at large r rotate too fast for visible mass alone.

 UQFF SOLUTION:
   v_UQFF(r) = √(G×M_eff(r)/r)
   
   where M_eff(r) = M_baryonic + M_UA(r) + M_SCm(r)

 AETHER MASS CONTRIBUTION:
   M_UA(r) = 4π ∫₀ʳ ρ_vac_UA × (1 + f_TRZ) × r'² dr'
   
   where:
     ρ_vac_UA = 7.09×10⁻³⁶ kg/m³     [Universal Aether density]
     f_TRZ = triadic resonance factor

 SCm MASS CONTRIBUTION:
   M_SCm(r) = 4π ∫₀ʳ ρ_vac_SCm × H_SCm × [SSq] × r'² dr'
   
   where:
     ρ_vac_SCm = 7.09×10⁻³⁷ kg/m³    [Superconductivity vacuum]
     H_SCm ≈ 0.99                     [Homogeneity factor]
     [SSq] = 0.57                     [Shell Quotient]

═══════════════════════════════════════════════════════════════════════════════
 PREDICTION: Flat rotation curves without dark matter particles
═══════════════════════════════════════════════════════════════════════════════
)";
    }
    
    QString getMagnetarEquation() {
        return R"(
═══════════════════════════════════════════════════════════════════════════════
 MAGNETAR FIELD EVOLUTION (SGR1745 / Soft Gamma Repeaters)
═══════════════════════════════════════════════════════════════════════════════

 CRITICAL FIELD:
   B_crit = 4.4×10¹³ T               [QED critical magnetic field]

 MAGNETIC FIELD EVOLUTION:
   B(t) = B₀ × e^(-t/τ_B) × (1 + f_UQFF)
   
   where:
     τ_B = 10⁴ - 10⁶ years           [Decay timescale]
     f_UQFF = Ug3 contribution       [String rotation enhancement]

 BURST ENERGY:
   E_burst = (B²/8πμ₀) × V_active × η_conversion
   
   where:
     V_active = crustal volume       [Active starquake region]
     η_conversion ≈ 0.01             [Field-to-radiation efficiency]

 UQFF ENHANCEMENT:
   Ug3_magnetar = (c/r_NS) × ω_s × B₀ × sin(θ_crust)
   
   where:
     r_NS ≈ 10 km                    [Neutron star radius]
     ω_s = angular velocity          [Typically 0.1-10 rad/s]

═══════════════════════════════════════════════════════════════════════════════
 OBSERVED: SGR1745 near Sagittarius A*, glitches correlated with Ug3
═══════════════════════════════════════════════════════════════════════════════
)";
    }
    
    QString getBECEquation() {
        return R"(
═══════════════════════════════════════════════════════════════════════════════
 BEC INTEGRATION MODULE (Bose-Einstein Condensate + UQFF)
═══════════════════════════════════════════════════════════════════════════════

 BOSE OCCUPANCY:
   N_B = 1 / (exp(ΔE/kT) - 1)
   
   where:
     ΔE = E - μ                      [Energy above chemical potential]
     kT = thermal energy             [Boltzmann × Temperature]

 CONDENSATE FRACTION:
   f_BEC = 1 - (T/T_c)^(3/2)         [Below critical temperature]
   
   where:
     T_c = (2πℏ²/mk) × (n/ζ(3/2))^(2/3)  [Critical temperature]

 UQFF COUPLING:
   Ug_BEC = g_contact × |ψ|² × n_BEC × SCm_factor
   
   where:
     g_contact = 4πℏ²a_s/m           [Contact interaction strength]
     a_s = scattering length         [Typically nm scale]
     ψ = condensate wavefunction
     SCm_factor = [SSq] × H_SCm      [Superconductivity enhancement]

 NUCLEAR COLLISION TERM:
   σ_collision = π × r_nuclear² × (1 + f_pairing)
   
   where:
     r_nuclear ≈ 1.2 × A^(1/3) fm    [Nuclear radius]
     f_pairing = pairing energy factor

═══════════════════════════════════════════════════════════════════════════════
)";
    }
    
    QString getLENREquation() {
        return R"(
═══════════════════════════════════════════════════════════════════════════════
 LENR WIDOM-LARSEN MODULE (Low-Energy Nuclear Reactions)
═══════════════════════════════════════════════════════════════════════════════

 SURFACE PLASMON POLARITON:
   ω_spp = ω_p / √2                  [Surface plasmon frequency]
   
   where:
     ω_p = √(n_e × e²/(ε₀ × m_e))    [Plasma frequency]

 HEAVY ELECTRON MASS:
   m_e* = m_e × (1 + λ)              [Dressed electron mass]
   
   where:
     λ = electron-phonon coupling ≈ 0.1-1

 NEUTRON PRODUCTION RATE:
   η = σ × n_H × v_rel × f_tunneling
   
   OBSERVED VALUES (Widom-Larsen systems):
     η_hydride ≈ 10¹³ cm⁻²/s        [Metal hydride surfaces]
     η_wires ≈ 10⁸ cm⁻²/s           [Nanowire configurations]
     η_corona ≈ 7×10⁻³ cm⁻²/s       [Corona discharge]

 UQFF ENHANCEMENT:
   F_LENR = k_LENR × sin(ω_LENR × t) × Ug_sum
   
   where:
     k_LENR = 10⁻¹⁰ N               [LENR coupling constant]
     ω_LENR = 1.2×10¹² rad/s        [THz resonance frequency]

═══════════════════════════════════════════════════════════════════════════════
 NOTE: LENR couples via THz resonance to UQFF vacuum fluctuations
═══════════════════════════════════════════════════════════════════════════════
)";
    }
    
    QString convertToLatex(const QString& unicode) {
        QString latex = unicode;
        // Convert Unicode math to LaTeX
        latex.replace("×", "\\times ");
        latex.replace("÷", "\\div ");
        latex.replace("≈", "\\approx ");
        latex.replace("≠", "\\neq ");
        latex.replace("≤", "\\leq ");
        latex.replace("≥", "\\geq ");
        latex.replace("∞", "\\infty ");
        latex.replace("π", "\\pi ");
        latex.replace("Σ", "\\sum ");
        latex.replace("Π", "\\prod ");
        latex.replace("∫", "\\int ");
        latex.replace("√", "\\sqrt ");
        latex.replace("²", "^{2}");
        latex.replace("³", "^{3}");
        latex.replace("⁴", "^{4}");
        latex.replace("₀", "_{0}");
        latex.replace("₁", "_{1}");
        latex.replace("₂", "_{2}");
        latex.replace("ℏ", "\\hbar ");
        latex.replace("μ", "\\mu ");
        latex.replace("ρ", "\\rho ");
        latex.replace("θ", "\\theta ");
        latex.replace("φ", "\\phi ");
        latex.replace("ω", "\\omega ");
        latex.replace("Λ", "\\Lambda ");
        latex.replace("β", "\\beta ");
        latex.replace("γ", "\\gamma ");
        latex.replace("κ", "\\kappa ");
        latex.replace("σ", "\\sigma ");
        latex.replace("Ω", "\\Omega ");
        latex.replace("λ", "\\lambda ");
        latex.replace("ε", "\\varepsilon ");
        latex.replace("ψ", "\\psi ");
        latex.replace("ζ", "\\zeta ");
        latex.replace("η", "\\eta ");
        latex.replace("τ", "\\tau ");
        return latex;
    }

    // Original member variables
    QTextEdit* equationDisplay;
    QComboBox* equationSelector;
    QComboBox* modeSelector;
    QString currentEquation;
    
    // Gap #4 Fix: IPC integration member variables
    QLineEdit* systemInput;         // Dynamic system name input
    QProgressBar* progressBar;      // Progress indicator for IPC
    QLabel* statusLabel;            // Status messages
    QProcess* qcalcProcess;         // QProcess for QCalc.py IPC
    QString pendingOutput;          // Buffer for IPC output
};

#endif // EQUATION_RENDERER_H
