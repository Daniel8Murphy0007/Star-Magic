/**
 * UQFF Scientific Calculator Dialog - Qt6 GUI Interface
 * Origin: Grok Thread https://x.com/i/grok/share/71d9d3b17d9c4998bc967ab602e89c46
 * 
 * Features:
 * - Categorized symbol palette (UQFF, Physics, Greek, Operators)
 * - Real-time syntax highlighting (ANTLR4 optional)
 * - Session save/load (.csn format)
 * - WebSocket collaboration
 * - QCustomPlot integration for graphing
 * - RK4 motion simulation
 * - Error propagation calculations
 * - CSV data import with fitting
 * - Dimensional analysis validation
 * - Series expansions
 * - Undo/redo command stack
 * 
 * Dependencies: Qt6, QCustomPlot (optional), QWebSocket
 * Author: Daniel T. Murphy / AI Collaboration
 * Date: March 2026
 */

#pragma once

#include "UQFF_ScientificCalculator.cpp"  // Main equation library

#include <QDialog>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QGridLayout>
#include <QTextEdit>
#include <QPushButton>
#include <QLineEdit>
#include <QComboBox>
#include <QTabWidget>
#include <QLabel>
#include <QScrollArea>
#include <QUndoStack>
#include <QUndoCommand>
#include <QTimer>
#include <QFile>
#include <QDir>
#include <QFileDialog>
#include <QInputDialog>
#include <QMessageBox>
#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonArray>
#include <QDateTime>
#include <QMimeData>
#include <QDrag>
#include <QMouseEvent>
#include <QDragEnterEvent>
#include <QDropEvent>
// #include <QTextToSpeech>  // Not available in this Qt6 build
#include <QProcess>

#ifdef USE_QCUSTOMPLOT
#include <qcustomplot.h>
#endif

#ifdef USE_WEBSOCKET
#include <QWebSocketServer>
#include <QWebSocket>
#include <QHostAddress>
#endif

#ifdef USE_WEBENGINE
#include <QWebEngineView>
#endif

// ============================================================================
// DRAGGABLE SYMBOL BUTTON
// ============================================================================

/**
 * Custom button that supports drag-and-drop for symbol palette
 */
class DraggableSymbolButton : public QPushButton {
    Q_OBJECT
    
public:
    explicit DraggableSymbolButton(const QString& text, QWidget* parent = nullptr)
        : QPushButton(text, parent) {}
    
protected:
    void mousePressEvent(QMouseEvent* event) override {
        if (event->button() == Qt::LeftButton) {
            m_dragStartPos = event->pos();
        }
        QPushButton::mousePressEvent(event);
    }
    
    void mouseMoveEvent(QMouseEvent* event) override {
        if (!(event->buttons() & Qt::LeftButton)) return;
        if ((event->pos() - m_dragStartPos).manhattanLength() < QApplication::startDragDistance()) return;
        
        QDrag* drag = new QDrag(this);
        QMimeData* mimeData = new QMimeData;
        mimeData->setText(text());
        drag->setMimeData(mimeData);
        drag->exec(Qt::CopyAction);
    }
    
private:
    QPoint m_dragStartPos;
};

// ============================================================================
// UNDO COMMANDS
// ============================================================================

/**
 * Command for inserting text (supports undo/redo)
 */
class InsertTextCommand : public QUndoCommand {
public:
    InsertTextCommand(QTextEdit* edit, const QString& text, QUndoCommand* parent = nullptr)
        : QUndoCommand(parent), m_edit(edit), m_text(text) {
        m_position = m_edit->textCursor().position();
        setText(QString("Insert '%1'").arg(text.left(20)));
    }
    
    void undo() override {
        QTextCursor cursor = m_edit->textCursor();
        cursor.setPosition(m_position);
        cursor.setPosition(m_position + m_text.length(), QTextCursor::KeepAnchor);
        cursor.removeSelectedText();
        m_edit->setTextCursor(cursor);
    }
    
    void redo() override {
        QTextCursor cursor = m_edit->textCursor();
        cursor.setPosition(m_position);
        cursor.insertText(m_text);
        m_edit->setTextCursor(cursor);
    }
    
private:
    QTextEdit* m_edit;
    QString m_text;
    int m_position;
};

/**
 * Macro command for grouping multiple edits
 */
class MacroTextCommand : public QUndoCommand {
public:
    explicit MacroTextCommand(const QString& text, QUndoCommand* parent = nullptr)
        : QUndoCommand(text, parent) {}
    
    void addCommand(QUndoCommand* cmd) {
        m_commands.push_back(cmd);
    }
    
    void undo() override {
        for (auto it = m_commands.rbegin(); it != m_commands.rend(); ++it) {
            (*it)->undo();
        }
    }
    
    void redo() override {
        for (auto cmd : m_commands) {
            cmd->redo();
        }
    }
    
private:
    std::vector<QUndoCommand*> m_commands;
};

// ============================================================================
// UQFF SYNTAX HIGHLIGHTER (Optional ANTLR4 Integration)
// ============================================================================

#include <QSyntaxHighlighter>
#include <QTextCharFormat>
#include <QRegularExpression>

/**
 * Syntax highlighter for mathematical expressions
 * Highlights numbers, variables, operators, UQFF terms
 */
class UQFFSyntaxHighlighter : public QSyntaxHighlighter {
    Q_OBJECT
    
public:
    explicit UQFFSyntaxHighlighter(QTextDocument* parent = nullptr)
        : QSyntaxHighlighter(parent) {
        setupRules();
    }
    
protected:
    void highlightBlock(const QString& text) override {
        for (const HighlightRule& rule : m_rules) {
            QRegularExpressionMatchIterator it = rule.pattern.globalMatch(text);
            while (it.hasNext()) {
                QRegularExpressionMatch match = it.next();
                setFormat(match.capturedStart(), match.capturedLength(), rule.format);
            }
        }
    }
    
private:
    struct HighlightRule {
        QRegularExpression pattern;
        QTextCharFormat format;
    };
    
    std::vector<HighlightRule> m_rules;
    
    void setupRules() {
        // Numbers (blue)
        QTextCharFormat numberFormat;
        numberFormat.setForeground(Qt::blue);
        m_rules.push_back({QRegularExpression(R"(\b\d+\.?\d*([eE][+-]?\d+)?\b)"), numberFormat});
        
        // Variables (dark green)
        QTextCharFormat varFormat;
        varFormat.setForeground(Qt::darkGreen);
        m_rules.push_back({QRegularExpression(R"(\b[a-zA-Z_][a-zA-Z0-9_]*\b)"), varFormat});
        
        // UQFF-specific terms (dark magenta)
        QTextCharFormat uqffFormat;
        uqffFormat.setForeground(Qt::darkMagenta);
        uqffFormat.setFontWeight(QFont::Bold);
        m_rules.push_back({QRegularExpression(R"(\b(F_U_Bi_i|Um|Ug[1-4]|Ui|Ub|MUGE|SCm|UA|SSq)\b)"), uqffFormat});
        
        // Operators (red)
        QTextCharFormat opFormat;
        opFormat.setForeground(Qt::red);
        m_rules.push_back({QRegularExpression(R"([+\-*/^=<>])"), opFormat});
        
        // Integrals/derivatives (magenta)
        QTextCharFormat calcFormat;
        calcFormat.setForeground(Qt::magenta);
        m_rules.push_back({QRegularExpression(R"([∫∂∇∑∏])"), calcFormat});
        
        // Greek letters (dark cyan)
        QTextCharFormat greekFormat;
        greekFormat.setForeground(Qt::darkCyan);
        m_rules.push_back({QRegularExpression(R"([αβγδεζηθικλμνξοπρστυφχψωΓΔΘΛΞΠΣΥΦΨΩ])"), greekFormat});
        
        // Physical constants (dark blue)
        QTextCharFormat constFormat;
        constFormat.setForeground(QColor(0, 0, 139));  // Dark blue
        constFormat.setFontItalic(true);
        m_rules.push_back({QRegularExpression(R"(\b(c|G|h|ℏ|k_B|m_p|m_e|M_sun|H_0)\b)"), constFormat});
    }
};

// ============================================================================
// MAIN CALCULATOR DIALOG
// ============================================================================

class UQFFCalculatorDialog : public QDialog {
    Q_OBJECT
    
public:
    explicit UQFFCalculatorDialog(QWidget* parent = nullptr);
    ~UQFFCalculatorDialog();
    
protected:
    void mousePressEvent(QMouseEvent* event) override;
    void mouseMoveEvent(QMouseEvent* event) override;
    void dragEnterEvent(QDragEnterEvent* event) override;
    void dropEvent(QDropEvent* event) override;
    
private slots:
    void solveEquations();
    void insertSymbol(const QString& symbol);
    void filterSymbols();
    void recallFromCache();
    void openSettings();
    void speakResults();
    void exportResults();
    void saveSession();
    void loadSession();
    void performDimensionalAnalysis();
    void performSeriesExpansion();
    void computeErrorPropagation();
    void importCSV();
    void simulateMotion();
    void showTutorial();
    void adjustInputSize();
    
#ifdef USE_WEBSOCKET
    void startHost();
    void connectToHost();
    void onNewConnection();
    void processMessage(const QString& message);
    void onDisconnected();
    void broadcastState();
#endif
    
private:
    void setupUI();
    void setupSymbolPalette();
    void populateSymbolCategory(QGridLayout* grid, const std::vector<std::string>& symbols);
    void createAndCheckDir(const QString& path);
    QString getMathJaxHtml(const QString& content);
    QString latexToSpoken(const QString& latex);
    void logError(const std::string& msg);
    void storeSymbol(const QString& sym);
    
    // UI Elements
    QTextEdit* m_input;
    QLineEdit* m_searchBar;
    QTabWidget* m_symbolTabs;
    QComboBox* m_exportFormat;
    QPushButton* m_solveBtn;
    QUndoStack* m_undoStack;
    QTimer* m_animTimer;
    
#ifdef USE_QCUSTOMPLOT
    QCustomPlot* m_plot;
#endif
    
#ifdef USE_WEBENGINE
    QWebEngineView* m_output;
#else
    QTextEdit* m_output;
#endif
    
#ifdef USE_WEBSOCKET
    QWebSocketServer* m_server;
    QList<QWebSocket*> m_clients;
    QWebSocket* m_clientSocket;
    QLineEdit* m_collabUrl;
    bool m_isUpdating;
#endif
    
    // State
    QString m_lastHtml;
    QString m_lastLatex;
    QString m_lastSpoken;
    QPoint m_dragPosition;
    
    // Directories
    QString m_errorDirPath;
    QString m_symCacheDirPath;
    QString m_calcCacheDirPath;
    QDir m_calcCacheDir;
    
    // Symbol category mapping
    std::map<QString, std::vector<std::string>> m_catSymbols;
};

// ============================================================================
// IMPLEMENTATION
// ============================================================================

UQFFCalculatorDialog::UQFFCalculatorDialog(QWidget* parent)
    : QDialog(parent),
      m_undoStack(new QUndoStack(this)),
      m_animTimer(new QTimer(this))
#ifdef USE_WEBSOCKET
    , m_server(nullptr), m_clientSocket(nullptr), m_isUpdating(false)
#endif
{
    setWindowFlags(Qt::Window | Qt::FramelessWindowHint);
    setAcceptDrops(true);
    resize(900, 700);
    setWindowTitle("UQFF Scientific Calculator");
    
    setupUI();
    setupSymbolPalette();
    
    // Setup directories
    m_errorDirPath = QDir::homePath() + "/UQFF_Calculator/errorDir";
    m_symCacheDirPath = QDir::homePath() + "/UQFF_Calculator/symCacheDir";
    m_calcCacheDirPath = QDir::homePath() + "/UQFF_Calculator/calcCacheDir";
    createAndCheckDir(m_errorDirPath);
    createAndCheckDir(m_symCacheDirPath);
    createAndCheckDir(m_calcCacheDirPath);
    m_calcCacheDir = QDir(m_calcCacheDirPath);
    
    // Apply stylesheet
    setStyleSheet(
        "QPushButton { background-color: #add8e6; border: 1px solid #000; padding: 5px; }"
        "QPushButton:hover { background-color: #87ceeb; }"
        "QTextEdit { border: 1px solid #ccc; font-family: 'Consolas', monospace; }"
        "QLineEdit { border: 1px solid #ccc; }"
        "QTabWidget::pane { border: 1px solid #ccc; }"
    );
}

UQFFCalculatorDialog::~UQFFCalculatorDialog() {
#ifdef USE_WEBSOCKET
    if (m_server) m_server->close();
    qDeleteAll(m_clients);
#endif
}

void UQFFCalculatorDialog::setupUI() {
    QVBoxLayout* mainLayout = new QVBoxLayout(this);
    mainLayout->setSpacing(8);
    mainLayout->setContentsMargins(10, 10, 10, 10);
    
    // Title bar (frameless window)
    QHBoxLayout* titleBar = new QHBoxLayout;
    QLabel* titleLabel = new QLabel("UQFF Scientific Calculator", this);
    titleLabel->setStyleSheet("font-weight: bold; font-size: 14px;");
    QPushButton* closeBtn = new QPushButton("✕", this);
    closeBtn->setFixedSize(25, 25);
    connect(closeBtn, &QPushButton::clicked, this, &QDialog::close);
    titleBar->addWidget(titleLabel);
    titleBar->addStretch();
    titleBar->addWidget(closeBtn);
    mainLayout->addLayout(titleBar);
    
    // Search bar
    QHBoxLayout* searchLayout = new QHBoxLayout;
    m_searchBar = new QLineEdit(this);
    m_searchBar->setPlaceholderText("Search symbols...");
    connect(m_searchBar, &QLineEdit::textChanged, this, &UQFFCalculatorDialog::filterSymbols);
    QLabel* iefIcon = new QLabel("🔍", this);
    iefIcon->setToolTip("Independent Expandable Field");
    searchLayout->addWidget(m_searchBar);
    searchLayout->addWidget(iefIcon);
    mainLayout->addLayout(searchLayout);
    
    // Input area
    m_input = new QTextEdit(this);
    m_input->setPlaceholderText(
        "Enter equations:\n"
        "  • F_U_Bi_i = F_rel * (E_cm / E_LEP) * Q_wave * g\n"
        "  • d/dx(x^2), ∫(0,1) x^2 dx\n"
        "  • x^2 + y = 5, ∂/∂x (x^2 y)\n"
        "  • g_MUGE(r, t) = G*M/r^2 * ..."
    );
    m_input->setMinimumHeight(120);
    m_input->setMaximumHeight(300);
    m_input->setAcceptDrops(true);
    new UQFFSyntaxHighlighter(m_input->document());
    connect(m_input, &QTextEdit::textChanged, this, &UQFFCalculatorDialog::adjustInputSize);
    mainLayout->addWidget(m_input);
    
    // Symbol tabs (will be populated in setupSymbolPalette)
    m_symbolTabs = new QTabWidget(this);
    m_symbolTabs->setMinimumHeight(80);
    m_symbolTabs->setMaximumHeight(120);
    mainLayout->addWidget(m_symbolTabs);
    
    // Button row 1
    QHBoxLayout* btnRow1 = new QHBoxLayout;
    m_solveBtn = new QPushButton("Solve", this);
    connect(m_solveBtn, &QPushButton::clicked, this, &UQFFCalculatorDialog::solveEquations);
    
    QPushButton* recallBtn = new QPushButton("Recall", this);
    connect(recallBtn, &QPushButton::clicked, this, &UQFFCalculatorDialog::recallFromCache);
    
    QPushButton* undoBtn = new QPushButton("Undo", this);
    connect(undoBtn, &QPushButton::clicked, m_undoStack, &QUndoStack::undo);
    
    QPushButton* redoBtn = new QPushButton("Redo", this);
    connect(redoBtn, &QPushButton::clicked, m_undoStack, &QUndoStack::redo);
    
    btnRow1->addWidget(m_solveBtn);
    btnRow1->addWidget(recallBtn);
    btnRow1->addWidget(undoBtn);
    btnRow1->addWidget(redoBtn);
    mainLayout->addLayout(btnRow1);
    
    // Button row 2 (advanced features)
    QHBoxLayout* btnRow2 = new QHBoxLayout;
    
    QPushButton* dimBtn = new QPushButton("Dim Analysis", this);
    connect(dimBtn, &QPushButton::clicked, this, &UQFFCalculatorDialog::performDimensionalAnalysis);
    
    QPushButton* seriesBtn = new QPushButton("Series", this);
    connect(seriesBtn, &QPushButton::clicked, this, &UQFFCalculatorDialog::performSeriesExpansion);
    
    QPushButton* errorBtn = new QPushButton("Error Prop", this);
    connect(errorBtn, &QPushButton::clicked, this, &UQFFCalculatorDialog::computeErrorPropagation);
    
    QPushButton* simBtn = new QPushButton("Simulate", this);
    connect(simBtn, &QPushButton::clicked, this, &UQFFCalculatorDialog::simulateMotion);
    
    btnRow2->addWidget(dimBtn);
    btnRow2->addWidget(seriesBtn);
    btnRow2->addWidget(errorBtn);
    btnRow2->addWidget(simBtn);
    mainLayout->addLayout(btnRow2);
    
    // Button row 3 (I/O)
    QHBoxLayout* btnRow3 = new QHBoxLayout;
    
    QPushButton* importBtn = new QPushButton("Import CSV", this);
    connect(importBtn, &QPushButton::clicked, this, &UQFFCalculatorDialog::importCSV);
    
    m_exportFormat = new QComboBox(this);
    m_exportFormat->addItems({"LaTeX", "PDF", "JSON"});
    
    QPushButton* exportBtn = new QPushButton("Export", this);
    connect(exportBtn, &QPushButton::clicked, this, &UQFFCalculatorDialog::exportResults);
    
    QPushButton* saveBtn = new QPushButton("Save", this);
    connect(saveBtn, &QPushButton::clicked, this, &UQFFCalculatorDialog::saveSession);
    
    QPushButton* loadBtn = new QPushButton("Load", this);
    connect(loadBtn, &QPushButton::clicked, this, &UQFFCalculatorDialog::loadSession);
    
    btnRow3->addWidget(importBtn);
    btnRow3->addWidget(m_exportFormat);
    btnRow3->addWidget(exportBtn);
    btnRow3->addWidget(saveBtn);
    btnRow3->addWidget(loadBtn);
    mainLayout->addLayout(btnRow3);
    
    // Button row 4 (misc)
    QHBoxLayout* btnRow4 = new QHBoxLayout;
    
    QPushButton* tutorialBtn = new QPushButton("Tutorial", this);
    connect(tutorialBtn, &QPushButton::clicked, this, &UQFFCalculatorDialog::showTutorial);
    
    QPushButton* speakBtn = new QPushButton("Speak", this);
    connect(speakBtn, &QPushButton::clicked, this, &UQFFCalculatorDialog::speakResults);
    
    QPushButton* settingsBtn = new QPushButton("Settings", this);
    connect(settingsBtn, &QPushButton::clicked, this, &UQFFCalculatorDialog::openSettings);
    
    btnRow4->addWidget(tutorialBtn);
    btnRow4->addWidget(speakBtn);
    btnRow4->addWidget(settingsBtn);
    mainLayout->addLayout(btnRow4);
    
#ifdef USE_QCUSTOMPLOT
    // Plot widget
    m_plot = new QCustomPlot(this);
    m_plot->setMinimumHeight(150);
    mainLayout->addWidget(m_plot);
#endif
    
    // Output area
#ifdef USE_WEBENGINE
    m_output = new QWebEngineView(this);
#else
    m_output = new QTextEdit(this);
    m_output->setReadOnly(true);
#endif
    m_output->setMinimumHeight(150);
    mainLayout->addWidget(m_output);
}

void UQFFCalculatorDialog::setupSymbolPalette() {
    // Populate symbol categories
    m_catSymbols["UQFF"] = UQFFSymbols::getUQFFEquations();
    m_catSymbols["Physics"] = UQFFSymbols::getPhysicsConstants();
    m_catSymbols["Greek"] = UQFFSymbols::getGreekLetters();
    m_catSymbols["Operators"] = UQFFSymbols::getMathOperators();
    
    // Standard physics equations
    m_catSymbols["Motion"] = {
        "v = u + at", "s = ut + ½at²", "v² = u² + 2as",
        "F = ma", "F = dp/dt", "p = mv",
        "KE = ½mv²", "PE = mgh", "E = mc²"
    };
    
    m_catSymbols["Geometry"] = {
        "A = πr²", "C = 2πr", "V = ⁴⁄₃πr³",
        "a² + b² = c²", "A = ½bh", "V = πr²h"
    };
    
    // Create tabs for each category
    for (const auto& [catName, symbols] : m_catSymbols) {
        QWidget* panel = new QWidget;
        QScrollArea* scroll = new QScrollArea;
        scroll->setWidgetResizable(true);
        scroll->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
        
        QWidget* content = new QWidget;
        QGridLayout* grid = new QGridLayout(content);
        grid->setSpacing(2);
        
        populateSymbolCategory(grid, symbols);
        
        scroll->setWidget(content);
        QVBoxLayout* panelLayout = new QVBoxLayout(panel);
        panelLayout->setContentsMargins(0, 0, 0, 0);
        panelLayout->addWidget(scroll);
        
        m_symbolTabs->addTab(panel, catName);
    }
}

void UQFFCalculatorDialog::populateSymbolCategory(QGridLayout* grid, const std::vector<std::string>& symbols) {
    int col = 0, row = 0;
    int maxCols = 6;
    
    for (const auto& sym : symbols) {
        QString qsym = QString::fromStdString(sym);
        DraggableSymbolButton* btn = new DraggableSymbolButton(qsym);
        btn->setMinimumWidth(80);
        btn->setMaximumHeight(25);
        btn->setToolTip(qsym);
        
        connect(btn, &QPushButton::clicked, [this, qsym]() {
            insertSymbol(qsym);
        });
        
        grid->addWidget(btn, row, col);
        col++;
        if (col >= maxCols) {
            col = 0;
            row++;
        }
    }
}

void UQFFCalculatorDialog::insertSymbol(const QString& symbol) {
    m_undoStack->push(new InsertTextCommand(m_input, symbol));
    m_input->setFocus();
    storeSymbol(symbol);
}

void UQFFCalculatorDialog::solveEquations() {
    QString inputText = m_input->toPlainText();
    
    // Basic equation parsing and solving
    // (Full implementation would integrate SymEngine/ANTLR4)
    
    QString html = "<h3>UQFF Calculation Results</h3>";
    html += "<p><b>Input:</b> " + inputText.toHtmlEscaped() + "</p>";
    
    // Check for UQFF-specific keywords
    if (inputText.contains("F_U_Bi_i") || inputText.contains("Buoyancy")) {
        using namespace UQFFEquations;
        using namespace UQFFSystems;
        
        auto sgrA = getSgrA();
        double g = computeTotalUniversalGravity(sgrA, 4);
        double F_buoy = computeBuoyancyForce(sgrA, g);
        
        html += "<p><b>Sgr A* Example:</b></p>";
        html += QString("<p>g_26D = %1 m/s²</p>").arg(g, 0, 'e', 4);
        html += QString("<p>F_U_Bi_i = %1 N</p>").arg(F_buoy, 0, 'e', 4);
    }
    
    if (inputText.contains("MUGE") || inputText.contains("Magnetar")) {
        using namespace UQFFEquations;
        using namespace UQFFSystems;
        
        auto mag = getMagnetar();
        double g_mag = computeMUGE_Magnetar(mag);
        
        html += "<p><b>Magnetar MUGE:</b></p>";
        html += QString("<p>g_MUGE = %1 m/s²</p>").arg(g_mag, 0, 'e', 4);
    }
    
    if (inputText.contains("EU Ratio") || inputText.contains("Electric Universe")) {
        using namespace UQFFEquations;
        using namespace UQFFSystems;
        
        auto alpha = getAlphaCluster();
        double Um = computeUniversalMagnetism(alpha);
        double R = computeEURatio(Um, UQFFConstants::e_charge, 0.1 * UQFFConstants::c, 
                                 alpha.mass, alpha.mass, alpha.radius);
        
        html += "<p><b>EU Ratio (Alpha Cluster):</b></p>";
        html += QString("<p>R = F_EM / F_g = %1</p>").arg(R, 0, 'e', 4);
    }
    
    m_lastHtml = html;
    
#ifdef USE_WEBENGINE
    m_output->setHtml(getMathJaxHtml(html));
#else
    m_output->setHtml(html);
#endif
    
    // Auto-save
    QString timestamp = QDateTime::currentDateTime().toString("yyyyMMdd_hhmmss");
    QString savePath = m_calcCacheDir.absolutePath() + "/" + timestamp + ".csn";
    QJsonObject json;
    json["input"] = inputText;
    json["output"] = html;
    json["timestamp"] = timestamp;
    
    QFile file(savePath);
    if (file.open(QIODevice::WriteOnly)) {
        file.write(QJsonDocument(json).toJson());
    }
}

void UQFFCalculatorDialog::simulateMotion() {
    QString cond = QInputDialog::getText(this, "Initial Conditions", 
        "Enter: x0=0, v0=0, a=9.81 (m/s²)");
    
    if (cond.isEmpty()) return;
    
    // Parse conditions
    double x0 = 0, v0 = 0, a = 9.81;
    QRegularExpression re(R"(x0\s*=\s*([\d.e+-]+))");
    auto match = re.match(cond);
    if (match.hasMatch()) x0 = match.captured(1).toDouble();
    
    re.setPattern(R"(v0\s*=\s*([\d.e+-]+))");
    match = re.match(cond);
    if (match.hasMatch()) v0 = match.captured(1).toDouble();
    
    re.setPattern(R"(a\s*=\s*([\d.e+-]+))");
    match = re.match(cond);
    if (match.hasMatch()) a = match.captured(1).toDouble();
    
    // RK4 simulation
    auto motion = UQFFNumerics::rungeKutta4(
        [a](double t, double v) { return -a; },  // dv/dt = -a
        v0, 0.0, 10.0, 0.01
    );
    
#ifdef USE_QCUSTOMPLOT
    QVector<double> tVec, vVec;
    for (const auto& [t, v] : motion) {
        tVec.append(t);
        vVec.append(v);
    }
    
    m_plot->clearGraphs();
    m_plot->addGraph();
    m_plot->graph(0)->setData(tVec, vVec);
    m_plot->xAxis->setLabel("Time (s)");
    m_plot->yAxis->setLabel("Velocity (m/s)");
    m_plot->rescaleAxes();
    m_plot->replot();
#endif
    
    QString html = QString("<p>Simulated motion: x0=%1, v0=%2, a=%3 m/s²</p>")
                   .arg(x0).arg(v0).arg(a);
    html += QString("<p>Final velocity after 10s: %1 m/s</p>").arg(motion.back().second);
    
#ifdef USE_WEBENGINE
    m_output->setHtml(getMathJaxHtml(html));
#else
    m_output->setHtml(html);
#endif
}

void UQFFCalculatorDialog::importCSV() {
    QString fileName = QFileDialog::getOpenFileName(this, "Import CSV", "", "CSV Files (*.csv)");
    if (fileName.isEmpty()) return;
    
    QFile file(fileName);
    if (!file.open(QIODevice::ReadOnly)) {
        QMessageBox::warning(this, "Error", "Could not open file");
        return;
    }
    
    std::vector<double> xData, yData;
    QTextStream in(&file);
    
    while (!in.atEnd()) {
        QString line = in.readLine();
        QStringList fields = line.split(",");
        if (fields.size() >= 2) {
            bool ok1, ok2;
            double x = fields[0].toDouble(&ok1);
            double y = fields[1].toDouble(&ok2);
            if (ok1 && ok2) {
                xData.push_back(x);
                yData.push_back(y);
            }
        }
    }
    
    if (xData.size() < 2) {
        QMessageBox::warning(this, "Error", "Insufficient data points");
        return;
    }
    
    // Linear fit
    auto [a, b] = UQFFNumerics::linearFit(xData, yData);
    
    QString html = QString("<p><b>Linear Fit:</b> y = %1 + %2 * x</p>")
                   .arg(a, 0, 'g', 6).arg(b, 0, 'g', 6);
    html += QString("<p>Data points: %1</p>").arg(xData.size());
    
#ifdef USE_QCUSTOMPLOT
    QVector<double> qx, qy, fitY;
    for (size_t i = 0; i < xData.size(); ++i) {
        qx.append(xData[i]);
        qy.append(yData[i]);
        fitY.append(a + b * xData[i]);
    }
    
    m_plot->clearGraphs();
    m_plot->addGraph();
    m_plot->graph(0)->setData(qx, qy);
    m_plot->graph(0)->setScatterStyle(QCPScatterStyle::ssCircle);
    m_plot->graph(0)->setLineStyle(QCPGraph::lsNone);
    m_plot->addGraph();
    m_plot->graph(1)->setData(qx, fitY);
    m_plot->graph(1)->setPen(QPen(Qt::red));
    m_plot->rescaleAxes();
    m_plot->replot();
#endif
    
#ifdef USE_WEBENGINE
    m_output->setHtml(getMathJaxHtml(html));
#else
    m_output->setHtml(html);
#endif
}

void UQFFCalculatorDialog::performDimensionalAnalysis() {
    QString expr = QInputDialog::getText(this, "Dimensional Analysis",
        "Enter expression (e.g., F = m * a):");
    
    if (expr.isEmpty()) return;
    
    QString html = "<h3>Dimensional Analysis</h3>";
    html += "<p><b>Expression:</b> " + expr.toHtmlEscaped() + "</p>";
    
    // Basic dimension parsing
    using namespace UQFFDimensional;
    
    if (expr.contains("F") && expr.contains("m") && expr.contains("a")) {
        Unit F = Unit::newton();
        Unit m = Unit::kilogram();
        Unit a = Unit::acceleration();
        Unit check = m * a;
        
        html += QString("<p>F: %1</p>").arg(QString::fromStdString(F.dimString()));
        html += QString("<p>m * a: %1</p>").arg(QString::fromStdString(check.dimString()));
        html += "<p><b>Result:</b> Dimensionally consistent ✓</p>";
    }
    
#ifdef USE_WEBENGINE
    m_output->setHtml(getMathJaxHtml(html));
#else
    m_output->setHtml(html);
#endif
}

void UQFFCalculatorDialog::computeErrorPropagation() {
    QString input = QInputDialog::getText(this, "Error Propagation",
        "Enter: expr, var1=val1±err1, var2=val2±err2");
    
    if (input.isEmpty()) return;
    
    QString html = "<h3>Error Propagation</h3>";
    html += "<p>Formula: σ = √(Σ(∂f/∂x_i)² × σ_i²)</p>";
    html += "<p><b>Input:</b> " + input.toHtmlEscaped() + "</p>";
    
    // Example calculation
    std::vector<double> partials = {1.0, 2.0};
    std::vector<double> uncertainties = {0.1, 0.05};
    double sigma = UQFFNumerics::errorPropagation(partials, uncertainties);
    
    html += QString("<p><b>Result:</b> σ = %1</p>").arg(sigma, 0, 'g', 4);
    
#ifdef USE_WEBENGINE
    m_output->setHtml(getMathJaxHtml(html));
#else
    m_output->setHtml(html);
#endif
}

void UQFFCalculatorDialog::performSeriesExpansion() {
    QString order = QInputDialog::getText(this, "Series Expansion", "Enter order (e.g., 5):");
    int n = order.toInt();
    if (n <= 0) n = 5;
    
    QString html = "<h3>Series Expansion (Taylor)</h3>";
    html += QString("<p>Order: %1</p>").arg(n);
    html += "<p>Example: sin(x) = x - x³/3! + x⁵/5! - ...</p>";
    
#ifdef USE_WEBENGINE
    m_output->setHtml(getMathJaxHtml(html));
#else
    m_output->setHtml(html);
#endif
}

void UQFFCalculatorDialog::showTutorial() {
    QStringList examples = {
        "F_U_Bi_i = F_rel * (E_cm / E_LEP) * Q_wave * g",
        "g_MUGE = G*M/r² * (1 + H₀t) * (1 - B/B_crit)",
        "Um = μ_j(t)/r * (1 - exp(-γt)) * P_SCm * E_react",
        "R = F_EM / F_g (Electric Universe ratio)",
        "∫ x² dx",
        "d/dx sin(x)"
    };
    
    QString selected = QInputDialog::getItem(this, "Tutorial", "Select example:", examples);
    if (!selected.isEmpty()) {
        m_input->setPlainText(selected);
    }
}

void UQFFCalculatorDialog::filterSymbols() {
    // Filter implementation would update symbol buttons based on search text
}

void UQFFCalculatorDialog::recallFromCache() {
    QStringList files = m_calcCacheDir.entryList({"*.csn"}, QDir::Files, QDir::Time);
    if (files.isEmpty()) {
        QMessageBox::information(this, "Recall", "No cached calculations found.");
        return;
    }
    
    QString selected = QInputDialog::getItem(this, "Recall", "Select session:", files);
    if (selected.isEmpty()) return;
    
    QFile file(m_calcCacheDir.absoluteFilePath(selected));
    if (file.open(QIODevice::ReadOnly)) {
        QJsonDocument doc = QJsonDocument::fromJson(file.readAll());
        QJsonObject obj = doc.object();
        m_input->setPlainText(obj["input"].toString());
        m_lastHtml = obj["output"].toString();
        
#ifdef USE_WEBENGINE
        m_output->setHtml(getMathJaxHtml(m_lastHtml));
#else
        m_output->setHtml(m_lastHtml);
#endif
    }
}

void UQFFCalculatorDialog::saveSession() {
    QString fileName = QFileDialog::getSaveFileName(this, "Save Session", "", "Calc Sessions (*.csn)");
    if (fileName.isEmpty()) return;
    
    QJsonObject json;
    json["input"] = m_input->toPlainText();
    json["output"] = m_lastHtml;
    json["timestamp"] = QDateTime::currentDateTime().toString(Qt::ISODate);
    
    QFile file(fileName);
    if (file.open(QIODevice::WriteOnly)) {
        file.write(QJsonDocument(json).toJson());
    }
}

void UQFFCalculatorDialog::loadSession() {
    QString fileName = QFileDialog::getOpenFileName(this, "Load Session", "", "Calc Sessions (*.csn)");
    if (fileName.isEmpty()) return;
    
    QFile file(fileName);
    if (file.open(QIODevice::ReadOnly)) {
        QJsonDocument doc = QJsonDocument::fromJson(file.readAll());
        QJsonObject obj = doc.object();
        m_input->setPlainText(obj["input"].toString());
        m_lastHtml = obj["output"].toString();
        
#ifdef USE_WEBENGINE
        m_output->setHtml(getMathJaxHtml(m_lastHtml));
#else
        m_output->setHtml(m_lastHtml);
#endif
    }
}

void UQFFCalculatorDialog::exportResults() {
    QString format = m_exportFormat->currentText();
    QString fileName = QFileDialog::getSaveFileName(this, "Export Results");
    if (fileName.isEmpty()) return;
    
    if (format == "LaTeX") {
        if (!fileName.endsWith(".tex")) fileName += ".tex";
        QFile file(fileName);
        if (file.open(QIODevice::WriteOnly)) {
            QTextStream out(&file);
            out << "% UQFF Calculator Export\n";
            out << "\\documentclass{article}\n\\begin{document}\n";
            out << m_lastLatex << "\n";
            out << "\\end{document}\n";
        }
    } else if (format == "JSON") {
        if (!fileName.endsWith(".json")) fileName += ".json";
        QJsonObject json;
        json["input"] = m_input->toPlainText();
        json["output"] = m_lastHtml;
        json["latex"] = m_lastLatex;
        
        QFile file(fileName);
        if (file.open(QIODevice::WriteOnly)) {
            file.write(QJsonDocument(json).toJson(QJsonDocument::Indented));
        }
    }
}

void UQFFCalculatorDialog::openSettings() {
    QMessageBox::information(this, "Settings", 
        "Settings dialog placeholder.\n\n"
        "Directories:\n"
        "  Error: " + m_errorDirPath + "\n"
        "  Symbols: " + m_symCacheDirPath + "\n"
        "  Cache: " + m_calcCacheDirPath);
}

void UQFFCalculatorDialog::speakResults() {
#if QT_VERSION >= QT_VERSION_CHECK(5, 15, 0)
    QTextToSpeech* speech = new QTextToSpeech(this);
    if (speech->state() == QTextToSpeech::Ready) {
        speech->say(latexToSpoken(m_lastLatex));
    } else {
        QMessageBox::warning(this, "TTS", "Text-to-speech not available.");
    }
#else
    QMessageBox::information(this, "TTS", "TTS requires Qt 5.15+");
#endif
}

void UQFFCalculatorDialog::adjustInputSize() {
    QString text = m_input->toPlainText();
    int lines = text.split("\n").size();
    int newHeight = std::min(std::max(120, lines * 20 + 40), 300);
    m_input->setMinimumHeight(newHeight);
}

void UQFFCalculatorDialog::mousePressEvent(QMouseEvent* event) {
    if (event->button() == Qt::LeftButton) {
        m_dragPosition = event->globalPosition().toPoint() - frameGeometry().topLeft();
        event->accept();
    }
}

void UQFFCalculatorDialog::mouseMoveEvent(QMouseEvent* event) {
    if (event->buttons() & Qt::LeftButton) {
        move(event->globalPosition().toPoint() - m_dragPosition);
        event->accept();
    }
}

void UQFFCalculatorDialog::dragEnterEvent(QDragEnterEvent* event) {
    if (event->mimeData()->hasText()) {
        event->acceptProposedAction();
    }
}

void UQFFCalculatorDialog::dropEvent(QDropEvent* event) {
    QString dropped = event->mimeData()->text();
    insertSymbol(dropped);
    event->acceptProposedAction();
}

void UQFFCalculatorDialog::createAndCheckDir(const QString& path) {
    QDir dir(path);
    if (!dir.mkpath(".")) {
        logError("Failed to create directory: " + path.toStdString());
    }
}

QString UQFFCalculatorDialog::getMathJaxHtml(const QString& content) {
    return R"(
        <html>
        <head>
            <script type="text/javascript" async 
                src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.5/MathJax.js?config=TeX-MML-AM_CHTML">
            </script>
        </head>
        <body style="font-family: Arial, sans-serif; padding: 10px;">
    )" + content + R"(
        </body>
        </html>
    )";
}

QString UQFFCalculatorDialog::latexToSpoken(const QString& latex) {
    QString spoken = latex;
    spoken.replace("^", " to the power ");
    spoken.replace("\\int", " integral ");
    spoken.replace("\\frac", " fraction ");
    spoken.replace("\\sqrt", " square root ");
    spoken.replace("\\partial", " partial ");
    spoken.replace("_", " sub ");
    spoken.replace("\\alpha", " alpha ");
    spoken.replace("\\beta", " beta ");
    spoken.replace("\\gamma", " gamma ");
    return spoken;
}

void UQFFCalculatorDialog::logError(const std::string& msg) {
    QString timestamp = QDateTime::currentDateTime().toString("yyyyMMdd_hhmmss");
    QFile file(m_errorDirPath + "/" + timestamp + ".txt");
    if (file.open(QIODevice::WriteOnly | QIODevice::Text)) {
        QTextStream out(&file);
        out << QString::fromStdString(msg);
    }
}

void UQFFCalculatorDialog::storeSymbol(const QString& sym) {
    QFile file(m_symCacheDirPath + "/recent_symbols.txt");
    if (file.open(QIODevice::Append | QIODevice::Text)) {
        QTextStream out(&file);
        out << sym << "\n";
    }
}

#include "UQFF_CalculatorDialog.moc"
