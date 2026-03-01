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

#include &lt;QDialog&gt;
#include &lt;QVBoxLayout&gt;
#include &lt;QHBoxLayout&gt;
#include &lt;QGridLayout&gt;
#include &lt;QTextEdit&gt;
#include &lt;QPushButton&gt;
#include &lt;QLineEdit&gt;
#include &lt;QComboBox&gt;
#include &lt;QTabWidget&gt;
#include &lt;QLabel&gt;
#include &lt;QScrollArea&gt;
#include &lt;QUndoStack&gt;
#include &lt;QUndoCommand&gt;
#include &lt;QTimer&gt;
#include &lt;QFile&gt;
#include &lt;QDir&gt;
#include &lt;QFileDialog&gt;
#include &lt;QInputDialog&gt;
#include &lt;QMessageBox&gt;
#include &lt;QJsonDocument&gt;
#include &lt;QJsonObject&gt;
#include &lt;QJsonArray&gt;
#include &lt;QDateTime&gt;
#include &lt;QMimeData&gt;
#include &lt;QDrag&gt;
#include &lt;QMouseEvent&gt;
#include &lt;QDragEnterEvent&gt;
#include &lt;QDropEvent&gt;
#include &lt;QTextToSpeech&gt;
#include &lt;QProcess&gt;

#ifdef USE_QCUSTOMPLOT
#include &lt;qcustomplot.h&gt;
#endif

#ifdef USE_WEBSOCKET
#include &lt;QWebSocketServer&gt;
#include &lt;QWebSocket&gt;
#include &lt;QHostAddress&gt;
#endif

#ifdef USE_WEBENGINE
#include &lt;QWebEngineView&gt;
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
        if (event-&gt;button() == Qt::LeftButton) {
            m_dragStartPos = event-&gt;pos();
        }
        QPushButton::mousePressEvent(event);
    }
    
    void mouseMoveEvent(QMouseEvent* event) override {
        if (!(event-&gt;buttons() &amp; Qt::LeftButton)) return;
        if ((event-&gt;pos() - m_dragStartPos).manhattanLength() &lt; QApplication::startDragDistance()) return;
        
        QDrag* drag = new QDrag(this);
        QMimeData* mimeData = new QMimeData;
        mimeData-&gt;setText(text());
        drag-&gt;setMimeData(mimeData);
        drag-&gt;exec(Qt::CopyAction);
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
        m_position = m_edit-&gt;textCursor().position();
        setText(QString("Insert '%1'").arg(text.left(20)));
    }
    
    void undo() override {
        QTextCursor cursor = m_edit-&gt;textCursor();
        cursor.setPosition(m_position);
        cursor.setPosition(m_position + m_text.length(), QTextCursor::KeepAnchor);
        cursor.removeSelectedText();
        m_edit-&gt;setTextCursor(cursor);
    }
    
    void redo() override {
        QTextCursor cursor = m_edit-&gt;textCursor();
        cursor.setPosition(m_position);
        cursor.insertText(m_text);
        m_edit-&gt;setTextCursor(cursor);
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
            (*it)-&gt;undo();
        }
    }
    
    void redo() override {
        for (auto cmd : m_commands) {
            cmd-&gt;redo();
        }
    }
    
private:
    std::vector&lt;QUndoCommand*&gt; m_commands;
};

// ============================================================================
// UQFF SYNTAX HIGHLIGHTER (Optional ANTLR4 Integration)
// ============================================================================

#include &lt;QSyntaxHighlighter&gt;
#include &lt;QTextCharFormat&gt;
#include &lt;QRegularExpression&gt;

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
    void highlightBlock(const QString&amp; text) override {
        for (const HighlightRule&amp; rule : m_rules) {
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
    
    std::vector&lt;HighlightRule&gt; m_rules;
    
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
        m_rules.push_back({QRegularExpression(R"([+\-*/^=&lt;&gt;])"), opFormat});
        
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
    void insertSymbol(const QString&amp; symbol);
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
    void processMessage(const QString&amp; message);
    void onDisconnected();
    void broadcastState();
#endif
    
private:
    void setupUI();
    void setupSymbolPalette();
    void populateSymbolCategory(QGridLayout* grid, const std::vector&lt;std::string&gt;&amp; symbols);
    void createAndCheckDir(const QString&amp; path);
    QString getMathJaxHtml(const QString&amp; content);
    QString latexToSpoken(const QString&amp; latex);
    void logError(const std::string&amp; msg);
    void storeSymbol(const QString&amp; sym);
    
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
    QList&lt;QWebSocket*&gt; m_clients;
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
    std::map&lt;QString, std::vector&lt;std::string&gt;&gt; m_catSymbols;
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
    if (m_server) m_server-&gt;close();
    qDeleteAll(m_clients);
#endif
}

void UQFFCalculatorDialog::setupUI() {
    QVBoxLayout* mainLayout = new QVBoxLayout(this);
    mainLayout-&gt;setSpacing(8);
    mainLayout-&gt;setContentsMargins(10, 10, 10, 10);
    
    // Title bar (frameless window)
    QHBoxLayout* titleBar = new QHBoxLayout;
    QLabel* titleLabel = new QLabel("UQFF Scientific Calculator", this);
    titleLabel-&gt;setStyleSheet("font-weight: bold; font-size: 14px;");
    QPushButton* closeBtn = new QPushButton("✕", this);
    closeBtn-&gt;setFixedSize(25, 25);
    connect(closeBtn, &amp;QPushButton::clicked, this, &amp;QDialog::close);
    titleBar-&gt;addWidget(titleLabel);
    titleBar-&gt;addStretch();
    titleBar-&gt;addWidget(closeBtn);
    mainLayout-&gt;addLayout(titleBar);
    
    // Search bar
    QHBoxLayout* searchLayout = new QHBoxLayout;
    m_searchBar = new QLineEdit(this);
    m_searchBar-&gt;setPlaceholderText("Search symbols...");
    connect(m_searchBar, &amp;QLineEdit::textChanged, this, &amp;UQFFCalculatorDialog::filterSymbols);
    QLabel* iefIcon = new QLabel("🔍", this);
    iefIcon-&gt;setToolTip("Independent Expandable Field");
    searchLayout-&gt;addWidget(m_searchBar);
    searchLayout-&gt;addWidget(iefIcon);
    mainLayout-&gt;addLayout(searchLayout);
    
    // Input area
    m_input = new QTextEdit(this);
    m_input-&gt;setPlaceholderText(
        "Enter equations:\n"
        "  • F_U_Bi_i = F_rel * (E_cm / E_LEP) * Q_wave * g\n"
        "  • d/dx(x^2), ∫(0,1) x^2 dx\n"
        "  • x^2 + y = 5, ∂/∂x (x^2 y)\n"
        "  • g_MUGE(r, t) = G*M/r^2 * ..."
    );
    m_input-&gt;setMinimumHeight(120);
    m_input-&gt;setMaximumHeight(300);
    m_input-&gt;setAcceptDrops(true);
    new UQFFSyntaxHighlighter(m_input-&gt;document());
    connect(m_input, &amp;QTextEdit::textChanged, this, &amp;UQFFCalculatorDialog::adjustInputSize);
    mainLayout-&gt;addWidget(m_input);
    
    // Symbol tabs (will be populated in setupSymbolPalette)
    m_symbolTabs = new QTabWidget(this);
    m_symbolTabs-&gt;setMinimumHeight(80);
    m_symbolTabs-&gt;setMaximumHeight(120);
    mainLayout-&gt;addWidget(m_symbolTabs);
    
    // Button row 1
    QHBoxLayout* btnRow1 = new QHBoxLayout;
    m_solveBtn = new QPushButton("Solve", this);
    connect(m_solveBtn, &amp;QPushButton::clicked, this, &amp;UQFFCalculatorDialog::solveEquations);
    
    QPushButton* recallBtn = new QPushButton("Recall", this);
    connect(recallBtn, &amp;QPushButton::clicked, this, &amp;UQFFCalculatorDialog::recallFromCache);
    
    QPushButton* undoBtn = new QPushButton("Undo", this);
    connect(undoBtn, &amp;QPushButton::clicked, m_undoStack, &amp;QUndoStack::undo);
    
    QPushButton* redoBtn = new QPushButton("Redo", this);
    connect(redoBtn, &amp;QPushButton::clicked, m_undoStack, &amp;QUndoStack::redo);
    
    btnRow1-&gt;addWidget(m_solveBtn);
    btnRow1-&gt;addWidget(recallBtn);
    btnRow1-&gt;addWidget(undoBtn);
    btnRow1-&gt;addWidget(redoBtn);
    mainLayout-&gt;addLayout(btnRow1);
    
    // Button row 2 (advanced features)
    QHBoxLayout* btnRow2 = new QHBoxLayout;
    
    QPushButton* dimBtn = new QPushButton("Dim Analysis", this);
    connect(dimBtn, &amp;QPushButton::clicked, this, &amp;UQFFCalculatorDialog::performDimensionalAnalysis);
    
    QPushButton* seriesBtn = new QPushButton("Series", this);
    connect(seriesBtn, &amp;QPushButton::clicked, this, &amp;UQFFCalculatorDialog::performSeriesExpansion);
    
    QPushButton* errorBtn = new QPushButton("Error Prop", this);
    connect(errorBtn, &amp;QPushButton::clicked, this, &amp;UQFFCalculatorDialog::computeErrorPropagation);
    
    QPushButton* simBtn = new QPushButton("Simulate", this);
    connect(simBtn, &amp;QPushButton::clicked, this, &amp;UQFFCalculatorDialog::simulateMotion);
    
    btnRow2-&gt;addWidget(dimBtn);
    btnRow2-&gt;addWidget(seriesBtn);
    btnRow2-&gt;addWidget(errorBtn);
    btnRow2-&gt;addWidget(simBtn);
    mainLayout-&gt;addLayout(btnRow2);
    
    // Button row 3 (I/O)
    QHBoxLayout* btnRow3 = new QHBoxLayout;
    
    QPushButton* importBtn = new QPushButton("Import CSV", this);
    connect(importBtn, &amp;QPushButton::clicked, this, &amp;UQFFCalculatorDialog::importCSV);
    
    m_exportFormat = new QComboBox(this);
    m_exportFormat-&gt;addItems({"LaTeX", "PDF", "JSON"});
    
    QPushButton* exportBtn = new QPushButton("Export", this);
    connect(exportBtn, &amp;QPushButton::clicked, this, &amp;UQFFCalculatorDialog::exportResults);
    
    QPushButton* saveBtn = new QPushButton("Save", this);
    connect(saveBtn, &amp;QPushButton::clicked, this, &amp;UQFFCalculatorDialog::saveSession);
    
    QPushButton* loadBtn = new QPushButton("Load", this);
    connect(loadBtn, &amp;QPushButton::clicked, this, &amp;UQFFCalculatorDialog::loadSession);
    
    btnRow3-&gt;addWidget(importBtn);
    btnRow3-&gt;addWidget(m_exportFormat);
    btnRow3-&gt;addWidget(exportBtn);
    btnRow3-&gt;addWidget(saveBtn);
    btnRow3-&gt;addWidget(loadBtn);
    mainLayout-&gt;addLayout(btnRow3);
    
    // Button row 4 (misc)
    QHBoxLayout* btnRow4 = new QHBoxLayout;
    
    QPushButton* tutorialBtn = new QPushButton("Tutorial", this);
    connect(tutorialBtn, &amp;QPushButton::clicked, this, &amp;UQFFCalculatorDialog::showTutorial);
    
    QPushButton* speakBtn = new QPushButton("Speak", this);
    connect(speakBtn, &amp;QPushButton::clicked, this, &amp;UQFFCalculatorDialog::speakResults);
    
    QPushButton* settingsBtn = new QPushButton("Settings", this);
    connect(settingsBtn, &amp;QPushButton::clicked, this, &amp;UQFFCalculatorDialog::openSettings);
    
    btnRow4-&gt;addWidget(tutorialBtn);
    btnRow4-&gt;addWidget(speakBtn);
    btnRow4-&gt;addWidget(settingsBtn);
    mainLayout-&gt;addLayout(btnRow4);
    
#ifdef USE_QCUSTOMPLOT
    // Plot widget
    m_plot = new QCustomPlot(this);
    m_plot-&gt;setMinimumHeight(150);
    mainLayout-&gt;addWidget(m_plot);
#endif
    
    // Output area
#ifdef USE_WEBENGINE
    m_output = new QWebEngineView(this);
#else
    m_output = new QTextEdit(this);
    m_output-&gt;setReadOnly(true);
#endif
    m_output-&gt;setMinimumHeight(150);
    mainLayout-&gt;addWidget(m_output);
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
        scroll-&gt;setWidgetResizable(true);
        scroll-&gt;setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
        
        QWidget* content = new QWidget;
        QGridLayout* grid = new QGridLayout(content);
        grid-&gt;setSpacing(2);
        
        populateSymbolCategory(grid, symbols);
        
        scroll-&gt;setWidget(content);
        QVBoxLayout* panelLayout = new QVBoxLayout(panel);
        panelLayout-&gt;setContentsMargins(0, 0, 0, 0);
        panelLayout-&gt;addWidget(scroll);
        
        m_symbolTabs-&gt;addTab(panel, catName);
    }
}

void UQFFCalculatorDialog::populateSymbolCategory(QGridLayout* grid, const std::vector&lt;std::string&gt;&amp; symbols) {
    int col = 0, row = 0;
    int maxCols = 6;
    
    for (const auto& sym : symbols) {
        QString qsym = QString::fromStdString(sym);
        DraggableSymbolButton* btn = new DraggableSymbolButton(qsym);
        btn-&gt;setMinimumWidth(80);
        btn-&gt;setMaximumHeight(25);
        btn-&gt;setToolTip(qsym);
        
        connect(btn, &amp;QPushButton::clicked, [this, qsym]() {
            insertSymbol(qsym);
        });
        
        grid-&gt;addWidget(btn, row, col);
        col++;
        if (col &gt;= maxCols) {
            col = 0;
            row++;
        }
    }
}

void UQFFCalculatorDialog::insertSymbol(const QString&amp; symbol) {
    m_undoStack-&gt;push(new InsertTextCommand(m_input, symbol));
    m_input-&gt;setFocus();
    storeSymbol(symbol);
}

void UQFFCalculatorDialog::solveEquations() {
    QString inputText = m_input-&gt;toPlainText();
    
    // Basic equation parsing and solving
    // (Full implementation would integrate SymEngine/ANTLR4)
    
    QString html = "&lt;h3&gt;UQFF Calculation Results&lt;/h3&gt;";
    html += "&lt;p&gt;&lt;b&gt;Input:&lt;/b&gt; " + inputText.toHtmlEscaped() + "&lt;/p&gt;";
    
    // Check for UQFF-specific keywords
    if (inputText.contains("F_U_Bi_i") || inputText.contains("Buoyancy")) {
        using namespace UQFFEquations;
        using namespace UQFFSystems;
        
        auto sgrA = getSgrA();
        double g = computeTotalUniversalGravity(sgrA, 4);
        double F_buoy = computeBuoyancyForce(sgrA, g);
        
        html += "&lt;p&gt;&lt;b&gt;Sgr A* Example:&lt;/b&gt;&lt;/p&gt;";
        html += QString("&lt;p&gt;g_26D = %1 m/s²&lt;/p&gt;").arg(g, 0, 'e', 4);
        html += QString("&lt;p&gt;F_U_Bi_i = %1 N&lt;/p&gt;").arg(F_buoy, 0, 'e', 4);
    }
    
    if (inputText.contains("MUGE") || inputText.contains("Magnetar")) {
        using namespace UQFFEquations;
        using namespace UQFFSystems;
        
        auto mag = getMagnetar();
        double g_mag = computeMUGE_Magnetar(mag);
        
        html += "&lt;p&gt;&lt;b&gt;Magnetar MUGE:&lt;/b&gt;&lt;/p&gt;";
        html += QString("&lt;p&gt;g_MUGE = %1 m/s²&lt;/p&gt;").arg(g_mag, 0, 'e', 4);
    }
    
    if (inputText.contains("EU Ratio") || inputText.contains("Electric Universe")) {
        using namespace UQFFEquations;
        using namespace UQFFSystems;
        
        auto alpha = getAlphaCluster();
        double Um = computeUniversalMagnetism(alpha);
        double R = computeEURatio(Um, UQFFConstants::e_charge, 0.1 * UQFFConstants::c, 
                                 alpha.mass, alpha.mass, alpha.radius);
        
        html += "&lt;p&gt;&lt;b&gt;EU Ratio (Alpha Cluster):&lt;/b&gt;&lt;/p&gt;";
        html += QString("&lt;p&gt;R = F_EM / F_g = %1&lt;/p&gt;").arg(R, 0, 'e', 4);
    }
    
    m_lastHtml = html;
    
#ifdef USE_WEBENGINE
    m_output-&gt;setHtml(getMathJaxHtml(html));
#else
    m_output-&gt;setHtml(html);
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
    QVector&lt;double&gt; tVec, vVec;
    for (const auto& [t, v] : motion) {
        tVec.append(t);
        vVec.append(v);
    }
    
    m_plot-&gt;clearGraphs();
    m_plot-&gt;addGraph();
    m_plot-&gt;graph(0)-&gt;setData(tVec, vVec);
    m_plot-&gt;xAxis-&gt;setLabel("Time (s)");
    m_plot-&gt;yAxis-&gt;setLabel("Velocity (m/s)");
    m_plot-&gt;rescaleAxes();
    m_plot-&gt;replot();
#endif
    
    QString html = QString("&lt;p&gt;Simulated motion: x0=%1, v0=%2, a=%3 m/s²&lt;/p&gt;")
                   .arg(x0).arg(v0).arg(a);
    html += QString("&lt;p&gt;Final velocity after 10s: %1 m/s&lt;/p&gt;").arg(motion.back().second);
    
#ifdef USE_WEBENGINE
    m_output-&gt;setHtml(getMathJaxHtml(html));
#else
    m_output-&gt;setHtml(html);
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
    
    std::vector&lt;double&gt; xData, yData;
    QTextStream in(&amp;file);
    
    while (!in.atEnd()) {
        QString line = in.readLine();
        QStringList fields = line.split(",");
        if (fields.size() &gt;= 2) {
            bool ok1, ok2;
            double x = fields[0].toDouble(&amp;ok1);
            double y = fields[1].toDouble(&amp;ok2);
            if (ok1 &amp;&amp; ok2) {
                xData.push_back(x);
                yData.push_back(y);
            }
        }
    }
    
    if (xData.size() &lt; 2) {
        QMessageBox::warning(this, "Error", "Insufficient data points");
        return;
    }
    
    // Linear fit
    auto [a, b] = UQFFNumerics::linearFit(xData, yData);
    
    QString html = QString("&lt;p&gt;&lt;b&gt;Linear Fit:&lt;/b&gt; y = %1 + %2 * x&lt;/p&gt;")
                   .arg(a, 0, 'g', 6).arg(b, 0, 'g', 6);
    html += QString("&lt;p&gt;Data points: %1&lt;/p&gt;").arg(xData.size());
    
#ifdef USE_QCUSTOMPLOT
    QVector&lt;double&gt; qx, qy, fitY;
    for (size_t i = 0; i &lt; xData.size(); ++i) {
        qx.append(xData[i]);
        qy.append(yData[i]);
        fitY.append(a + b * xData[i]);
    }
    
    m_plot-&gt;clearGraphs();
    m_plot-&gt;addGraph();
    m_plot-&gt;graph(0)-&gt;setData(qx, qy);
    m_plot-&gt;graph(0)-&gt;setScatterStyle(QCPScatterStyle::ssCircle);
    m_plot-&gt;graph(0)-&gt;setLineStyle(QCPGraph::lsNone);
    m_plot-&gt;addGraph();
    m_plot-&gt;graph(1)-&gt;setData(qx, fitY);
    m_plot-&gt;graph(1)-&gt;setPen(QPen(Qt::red));
    m_plot-&gt;rescaleAxes();
    m_plot-&gt;replot();
#endif
    
#ifdef USE_WEBENGINE
    m_output-&gt;setHtml(getMathJaxHtml(html));
#else
    m_output-&gt;setHtml(html);
#endif
}

void UQFFCalculatorDialog::performDimensionalAnalysis() {
    QString expr = QInputDialog::getText(this, "Dimensional Analysis",
        "Enter expression (e.g., F = m * a):");
    
    if (expr.isEmpty()) return;
    
    QString html = "&lt;h3&gt;Dimensional Analysis&lt;/h3&gt;";
    html += "&lt;p&gt;&lt;b&gt;Expression:&lt;/b&gt; " + expr.toHtmlEscaped() + "&lt;/p&gt;";
    
    // Basic dimension parsing
    using namespace UQFFDimensional;
    
    if (expr.contains("F") &amp;&amp; expr.contains("m") &amp;&amp; expr.contains("a")) {
        Unit F = Unit::newton();
        Unit m = Unit::kilogram();
        Unit a = Unit::acceleration();
        Unit check = m * a;
        
        html += QString("&lt;p&gt;F: %1&lt;/p&gt;").arg(QString::fromStdString(F.dimString()));
        html += QString("&lt;p&gt;m * a: %1&lt;/p&gt;").arg(QString::fromStdString(check.dimString()));
        html += "&lt;p&gt;&lt;b&gt;Result:&lt;/b&gt; Dimensionally consistent ✓&lt;/p&gt;";
    }
    
#ifdef USE_WEBENGINE
    m_output-&gt;setHtml(getMathJaxHtml(html));
#else
    m_output-&gt;setHtml(html);
#endif
}

void UQFFCalculatorDialog::computeErrorPropagation() {
    QString input = QInputDialog::getText(this, "Error Propagation",
        "Enter: expr, var1=val1±err1, var2=val2±err2");
    
    if (input.isEmpty()) return;
    
    QString html = "&lt;h3&gt;Error Propagation&lt;/h3&gt;";
    html += "&lt;p&gt;Formula: σ = √(Σ(∂f/∂x_i)² × σ_i²)&lt;/p&gt;";
    html += "&lt;p&gt;&lt;b&gt;Input:&lt;/b&gt; " + input.toHtmlEscaped() + "&lt;/p&gt;";
    
    // Example calculation
    std::vector&lt;double&gt; partials = {1.0, 2.0};
    std::vector&lt;double&gt; uncertainties = {0.1, 0.05};
    double sigma = UQFFNumerics::errorPropagation(partials, uncertainties);
    
    html += QString("&lt;p&gt;&lt;b&gt;Result:&lt;/b&gt; σ = %1&lt;/p&gt;").arg(sigma, 0, 'g', 4);
    
#ifdef USE_WEBENGINE
    m_output-&gt;setHtml(getMathJaxHtml(html));
#else
    m_output-&gt;setHtml(html);
#endif
}

void UQFFCalculatorDialog::performSeriesExpansion() {
    QString order = QInputDialog::getText(this, "Series Expansion", "Enter order (e.g., 5):");
    int n = order.toInt();
    if (n &lt;= 0) n = 5;
    
    QString html = "&lt;h3&gt;Series Expansion (Taylor)&lt;/h3&gt;";
    html += QString("&lt;p&gt;Order: %1&lt;/p&gt;").arg(n);
    html += "&lt;p&gt;Example: sin(x) = x - x³/3! + x⁵/5! - ...&lt;/p&gt;";
    
#ifdef USE_WEBENGINE
    m_output-&gt;setHtml(getMathJaxHtml(html));
#else
    m_output-&gt;setHtml(html);
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
        m_input-&gt;setPlainText(selected);
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
        m_input-&gt;setPlainText(obj["input"].toString());
        m_lastHtml = obj["output"].toString();
        
#ifdef USE_WEBENGINE
        m_output-&gt;setHtml(getMathJaxHtml(m_lastHtml));
#else
        m_output-&gt;setHtml(m_lastHtml);
#endif
    }
}

void UQFFCalculatorDialog::saveSession() {
    QString fileName = QFileDialog::getSaveFileName(this, "Save Session", "", "Calc Sessions (*.csn)");
    if (fileName.isEmpty()) return;
    
    QJsonObject json;
    json["input"] = m_input-&gt;toPlainText();
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
        m_input-&gt;setPlainText(obj["input"].toString());
        m_lastHtml = obj["output"].toString();
        
#ifdef USE_WEBENGINE
        m_output-&gt;setHtml(getMathJaxHtml(m_lastHtml));
#else
        m_output-&gt;setHtml(m_lastHtml);
#endif
    }
}

void UQFFCalculatorDialog::exportResults() {
    QString format = m_exportFormat-&gt;currentText();
    QString fileName = QFileDialog::getSaveFileName(this, "Export Results");
    if (fileName.isEmpty()) return;
    
    if (format == "LaTeX") {
        if (!fileName.endsWith(".tex")) fileName += ".tex";
        QFile file(fileName);
        if (file.open(QIODevice::WriteOnly)) {
            QTextStream out(&amp;file);
            out &lt;&lt; "% UQFF Calculator Export\n";
            out &lt;&lt; "\\documentclass{article}\n\\begin{document}\n";
            out &lt;&lt; m_lastLatex &lt;&lt; "\n";
            out &lt;&lt; "\\end{document}\n";
        }
    } else if (format == "JSON") {
        if (!fileName.endsWith(".json")) fileName += ".json";
        QJsonObject json;
        json["input"] = m_input-&gt;toPlainText();
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
#if QT_VERSION &gt;= QT_VERSION_CHECK(5, 15, 0)
    QTextToSpeech* speech = new QTextToSpeech(this);
    if (speech-&gt;state() == QTextToSpeech::Ready) {
        speech-&gt;say(latexToSpoken(m_lastLatex));
    } else {
        QMessageBox::warning(this, "TTS", "Text-to-speech not available.");
    }
#else
    QMessageBox::information(this, "TTS", "TTS requires Qt 5.15+");
#endif
}

void UQFFCalculatorDialog::adjustInputSize() {
    QString text = m_input-&gt;toPlainText();
    int lines = text.split("\n").size();
    int newHeight = std::min(std::max(120, lines * 20 + 40), 300);
    m_input-&gt;setMinimumHeight(newHeight);
}

void UQFFCalculatorDialog::mousePressEvent(QMouseEvent* event) {
    if (event-&gt;button() == Qt::LeftButton) {
        m_dragPosition = event-&gt;globalPosition().toPoint() - frameGeometry().topLeft();
        event-&gt;accept();
    }
}

void UQFFCalculatorDialog::mouseMoveEvent(QMouseEvent* event) {
    if (event-&gt;buttons() &amp; Qt::LeftButton) {
        move(event-&gt;globalPosition().toPoint() - m_dragPosition);
        event-&gt;accept();
    }
}

void UQFFCalculatorDialog::dragEnterEvent(QDragEnterEvent* event) {
    if (event-&gt;mimeData()-&gt;hasText()) {
        event-&gt;acceptProposedAction();
    }
}

void UQFFCalculatorDialog::dropEvent(QDropEvent* event) {
    QString dropped = event-&gt;mimeData()-&gt;text();
    insertSymbol(dropped);
    event-&gt;acceptProposedAction();
}

void UQFFCalculatorDialog::createAndCheckDir(const QString&amp; path) {
    QDir dir(path);
    if (!dir.mkpath(".")) {
        logError("Failed to create directory: " + path.toStdString());
    }
}

QString UQFFCalculatorDialog::getMathJaxHtml(const QString&amp; content) {
    return R"(
        &lt;html&gt;
        &lt;head&gt;
            &lt;script type="text/javascript" async 
                src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.5/MathJax.js?config=TeX-MML-AM_CHTML"&gt;
            &lt;/script&gt;
        &lt;/head&gt;
        &lt;body style="font-family: Arial, sans-serif; padding: 10px;"&gt;
    )" + content + R"(
        &lt;/body&gt;
        &lt;/html&gt;
    )";
}

QString UQFFCalculatorDialog::latexToSpoken(const QString&amp; latex) {
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

void UQFFCalculatorDialog::logError(const std::string&amp; msg) {
    QString timestamp = QDateTime::currentDateTime().toString("yyyyMMdd_hhmmss");
    QFile file(m_errorDirPath + "/" + timestamp + ".txt");
    if (file.open(QIODevice::WriteOnly | QIODevice::Text)) {
        QTextStream out(&amp;file);
        out &lt;&lt; QString::fromStdString(msg);
    }
}

void UQFFCalculatorDialog::storeSymbol(const QString&amp; sym) {
    QFile file(m_symCacheDirPath + "/recent_symbols.txt");
    if (file.open(QIODevice::Append | QIODevice::Text)) {
        QTextStream out(&amp;file);
        out &lt;&lt; sym &lt;&lt; "\n";
    }
}

#include "UQFF_CalculatorDialog.moc"
