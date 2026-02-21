// ═══════════════════════════════════════════════════════════════════════════════
// SOURCE2_WIDGETS_ENHANCED.H - Enhanced Widgets for Source2 GUI
// ═══════════════════════════════════════════════════════════════════════════════
// 
// Contains:
//   1. SessionLogWidget - Real-time log viewer (Tab 9)
//   2. PythonBridge - Bidirectional C++ ↔ Python communication
//   3. ComparisonDashboard - Split-view results comparison
//   4. SessionPersistence - Save/restore session state
//
// © 2025-2026 Daniel T. Murphy - Star-Magic UQFF Framework
// ═══════════════════════════════════════════════════════════════════════════════

#ifndef SOURCE2_WIDGETS_ENHANCED_H
#define SOURCE2_WIDGETS_ENHANCED_H

#include <QWidget>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QTextEdit>
#include <QLineEdit>
#include <QPushButton>
#include <QComboBox>
#include <QCheckBox>
#include <QLabel>
#include <QGroupBox>
#include <QSplitter>
#include <QTableWidget>
#include <QHeaderView>
#include <QProcess>
#include <QFile>
#include <QFileDialog>
#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonArray>
#include <QDateTime>
#include <QTimer>
#include <QScrollBar>
#include <QRegularExpression>

#include "source2_event_bus.h"

// ═══════════════════════════════════════════════════════════════════════════════
// 1. SESSION LOG WIDGET - Real-Time Log Viewer (Tab 9)
// ═══════════════════════════════════════════════════════════════════════════════

class SessionLogWidget : public QWidget {
    Q_OBJECT
    
public:
    SessionLogWidget(QWidget* parent = nullptr) : QWidget(parent) {
        setupUI();
        connectEventBus();
        
        // Log startup
        UQFF_LOG_INFO("SessionLog", "Session log viewer initialized");
    }
    
private:
    void setupUI() {
        QVBoxLayout* mainLayout = new QVBoxLayout(this);
        mainLayout->setSpacing(5);
        mainLayout->setContentsMargins(10, 10, 10, 10);
        
        // ═══════════════════════════════════════════════════════════════════════
        // HEADER
        // ═══════════════════════════════════════════════════════════════════════
        
        QLabel* header = new QLabel(
            "<h2 style='color: #4CAF50;'>📋 Session Log Viewer</h2>"
            "<p style='color: #888;'>Real-time logging from all Source2 components</p>"
        );
        mainLayout->addWidget(header);
        
        // ═══════════════════════════════════════════════════════════════════════
        // FILTER CONTROLS
        // ═══════════════════════════════════════════════════════════════════════
        
        QGroupBox* filterGroup = new QGroupBox("Filters");
        QHBoxLayout* filterLayout = new QHBoxLayout(filterGroup);
        
        // Level filter
        filterLayout->addWidget(new QLabel("Level:"));
        levelFilter = new QComboBox();
        levelFilter->addItem("All", -1);
        levelFilter->addItem("DEBUG", static_cast<int>(LogLevel::DEBUG));
        levelFilter->addItem("INFO", static_cast<int>(LogLevel::INFO));
        levelFilter->addItem("PHYSICS", static_cast<int>(LogLevel::PHYSICS));
        levelFilter->addItem("VALIDATION", static_cast<int>(LogLevel::VALIDATION));
        levelFilter->addItem("WARNING", static_cast<int>(LogLevel::WARNING));
        levelFilter->addItem("ERROR", static_cast<int>(LogLevel::ERROR));
        connect(levelFilter, QOverload<int>::of(&QComboBox::currentIndexChanged),
                this, &SessionLogWidget::applyFilters);
        filterLayout->addWidget(levelFilter);
        
        // Source filter
        filterLayout->addWidget(new QLabel("Source:"));
        sourceFilter = new QComboBox();
        sourceFilter->addItem("All");
        sourceFilter->addItem("C++");
        sourceFilter->addItem("Python");
        sourceFilter->addItem("Grok");
        sourceFilter->addItem("Ollama");
        sourceFilter->addItem("API");
        sourceFilter->addItem("Validation");
        sourceFilter->setEditable(true);
        connect(sourceFilter, &QComboBox::currentTextChanged,
                this, &SessionLogWidget::applyFilters);
        filterLayout->addWidget(sourceFilter);
        
        // Search box
        filterLayout->addWidget(new QLabel("Search:"));
        searchBox = new QLineEdit();
        searchBox->setPlaceholderText("Filter by text...");
        connect(searchBox, &QLineEdit::textChanged, this, &SessionLogWidget::applyFilters);
        filterLayout->addWidget(searchBox);
        
        // Auto-scroll checkbox
        autoScrollCheck = new QCheckBox("Auto-scroll");
        autoScrollCheck->setChecked(true);
        filterLayout->addWidget(autoScrollCheck);
        
        filterLayout->addStretch();
        mainLayout->addWidget(filterGroup);
        
        // ═══════════════════════════════════════════════════════════════════════
        // LOG DISPLAY
        // ═══════════════════════════════════════════════════════════════════════
        
        logDisplay = new QTextEdit();
        logDisplay->setReadOnly(true);
        logDisplay->setStyleSheet(
            "QTextEdit {"
            "  background-color: #1a1a1a;"
            "  color: #e0e0e0;"
            "  font-family: 'Consolas', 'Monaco', monospace;"
            "  font-size: 11px;"
            "  border: 1px solid #333;"
            "}"
        );
        mainLayout->addWidget(logDisplay, 1);
        
        // ═══════════════════════════════════════════════════════════════════════
        // STATISTICS BAR
        // ═══════════════════════════════════════════════════════════════════════
        
        QHBoxLayout* statsLayout = new QHBoxLayout();
        
        statsLabel = new QLabel("Entries: 0 | Physics: 0 | Validations: 0 | Errors: 0");
        statsLabel->setStyleSheet("color: #888;");
        statsLayout->addWidget(statsLabel);
        
        statsLayout->addStretch();
        
        // Action buttons
        QPushButton* clearBtn = new QPushButton("🗑️ Clear");
        connect(clearBtn, &QPushButton::clicked, this, &SessionLogWidget::clearLogs);
        statsLayout->addWidget(clearBtn);
        
        QPushButton* exportBtn = new QPushButton("📥 Export");
        connect(exportBtn, &QPushButton::clicked, this, &SessionLogWidget::exportLogs);
        statsLayout->addWidget(exportBtn);
        
        QPushButton* pauseBtn = new QPushButton("⏸️ Pause");
        pauseBtn->setCheckable(true);
        connect(pauseBtn, &QPushButton::toggled, this, &SessionLogWidget::togglePause);
        statsLayout->addWidget(pauseBtn);
        
        mainLayout->addLayout(statsLayout);
        
        setLayout(mainLayout);
    }
    
    void connectEventBus() {
        connect(&UQFFEventBus::instance(), &UQFFEventBus::logMessage,
                this, &SessionLogWidget::onLogMessage);
    }
    
private slots:
    void onLogMessage(const LogEntry& entry) {
        if (m_paused) return;
        
        m_allEntries.append(entry);
        updateStats();
        
        // Check filters
        if (!passesFilter(entry)) return;
        
        // Format and display
        QString html = formatEntry(entry);
        logDisplay->append(html);
        
        // Auto-scroll
        if (autoScrollCheck->isChecked()) {
            QScrollBar* sb = logDisplay->verticalScrollBar();
            sb->setValue(sb->maximum());
        }
    }
    
    void applyFilters() {
        logDisplay->clear();
        for (const auto& entry : m_allEntries) {
            if (passesFilter(entry)) {
                logDisplay->append(formatEntry(entry));
            }
        }
    }
    
    void clearLogs() {
        logDisplay->clear();
        m_allEntries.clear();
        UQFFEventBus::instance().clearLogHistory();
        updateStats();
        UQFF_LOG_INFO("SessionLog", "Log history cleared");
    }
    
    void exportLogs() {
        QString filename = QFileDialog::getSaveFileName(
            this, "Export Logs", 
            QString("uqff_session_%1.log").arg(QDateTime::currentDateTime().toString("yyyyMMdd_hhmmss")),
            "Log Files (*.log);;JSON Files (*.json);;All Files (*.*)"
        );
        
        if (filename.isEmpty()) return;
        
        QFile file(filename);
        if (file.open(QIODevice::WriteOnly | QIODevice::Text)) {
            QTextStream out(&file);
            
            if (filename.endsWith(".json")) {
                // JSON export
                QJsonArray arr;
                for (const auto& entry : m_allEntries) {
                    arr.append(entry.toJson());
                }
                out << QJsonDocument(arr).toJson(QJsonDocument::Indented);
            } else {
                // Plain text export
                for (const auto& entry : m_allEntries) {
                    out << QString("[%1] [%2] [%3] %4\n")
                           .arg(entry.timestamp.toString("yyyy-MM-dd hh:mm:ss.zzz"))
                           .arg(logLevelToString(entry.level), -10)
                           .arg(entry.source, -15)
                           .arg(entry.message);
                }
            }
            
            file.close();
            UQFF_LOG_INFO("SessionLog", QString("Logs exported to: %1").arg(filename));
        }
    }
    
    void togglePause(bool paused) {
        m_paused = paused;
    }
    
private:
    bool passesFilter(const LogEntry& entry) {
        // Level filter
        int levelValue = levelFilter->currentData().toInt();
        if (levelValue >= 0 && static_cast<int>(entry.level) != levelValue) {
            return false;
        }
        
        // Source filter
        QString sourceText = sourceFilter->currentText();
        if (sourceText != "All" && !entry.source.contains(sourceText, Qt::CaseInsensitive)) {
            return false;
        }
        
        // Search filter
        QString searchText = searchBox->text();
        if (!searchText.isEmpty() && 
            !entry.message.contains(searchText, Qt::CaseInsensitive) &&
            !entry.source.contains(searchText, Qt::CaseInsensitive)) {
            return false;
        }
        
        return true;
    }
    
    QString formatEntry(const LogEntry& entry) {
        QString color = logLevelToColor(entry.level);
        QString levelBadge = QString("<span style='background-color: %1; color: white; "
                                     "padding: 1px 4px; border-radius: 2px; font-size: 9px;'>%2</span>")
                             .arg(color).arg(logLevelToString(entry.level));
        
        return QString("<p style='margin: 2px 0;'>"
                      "<span style='color: #666;'>%1</span> "
                      "%2 "
                      "<span style='color: #2196F3;'>[%3]</span> "
                      "<span style='color: #e0e0e0;'>%4</span>"
                      "</p>")
               .arg(entry.timestamp.toString("hh:mm:ss.zzz"))
               .arg(levelBadge)
               .arg(entry.source)
               .arg(entry.message.toHtmlEscaped());
    }
    
    void updateStats() {
        int total = m_allEntries.size();
        int physics = 0, validations = 0, errors = 0;
        
        for (const auto& e : m_allEntries) {
            if (e.level == LogLevel::PHYSICS) physics++;
            else if (e.level == LogLevel::VALIDATION) validations++;
            else if (e.level == LogLevel::ERROR) errors++;
        }
        
        statsLabel->setText(QString("Entries: %1 | Physics: %2 | Validations: %3 | Errors: %4")
                           .arg(total).arg(physics).arg(validations).arg(errors));
    }
    
    QTextEdit* logDisplay;
    QComboBox* levelFilter;
    QComboBox* sourceFilter;
    QLineEdit* searchBox;
    QCheckBox* autoScrollCheck;
    QLabel* statsLabel;
    
    QList<LogEntry> m_allEntries;
    bool m_paused = false;
};


// ═══════════════════════════════════════════════════════════════════════════════
// 2. PYTHON BRIDGE - Bidirectional C++ ↔ Python Communication
// ═══════════════════════════════════════════════════════════════════════════════

class PythonBridge : public QObject {
    Q_OBJECT
    
public:
    PythonBridge(QObject* parent = nullptr) : QObject(parent) {
        m_process = new QProcess(this);
        
        connect(m_process, &QProcess::readyReadStandardOutput,
                this, &PythonBridge::onStdoutReady);
        connect(m_process, &QProcess::readyReadStandardError,
                this, &PythonBridge::onStderrReady);
        connect(m_process, QOverload<int, QProcess::ExitStatus>::of(&QProcess::finished),
                this, &PythonBridge::onProcessFinished);
    }
    
    ~PythonBridge() {
        stop();
    }
    
    // Start the Python JSON-RPC server
    bool start(const QString& scriptPath = "CondensedPhysics.py") {
        if (m_process->state() != QProcess::NotRunning) {
            return false;
        }
        
        QString python = findPythonExecutable();
        QStringList args;
        args << "-u" << scriptPath << "--json-rpc";
        
        m_process->start(python, args);
        bool started = m_process->waitForStarted(5000);
        
        if (started) {
            m_running = true;
            UQFF_LOG_INFO("PythonBridge", QString("Started Python bridge: %1").arg(scriptPath));
        } else {
            UQFF_LOG_ERROR("PythonBridge", QString("Failed to start Python: %1").arg(m_process->errorString()));
        }
        
        return started;
    }
    
    void stop() {
        if (m_process->state() != QProcess::NotRunning) {
            m_process->write("{\"method\": \"__shutdown__\"}\n");
            m_process->waitForFinished(3000);
            if (m_process->state() != QProcess::NotRunning) {
                m_process->kill();
            }
        }
        m_running = false;
    }
    
    bool isRunning() const { return m_running; }
    
    // ═══════════════════════════════════════════════════════════════════════════
    // REMOTE PROCEDURE CALLS
    // ═══════════════════════════════════════════════════════════════════════════
    
    void callAsync(const QString& method, const QJsonObject& params) {
        if (!m_running) {
            UQFF_LOG_WARNING("PythonBridge", "Cannot call method - bridge not running");
            return;
        }
        
        QJsonObject request;
        request["id"] = m_requestId++;
        request["method"] = method;
        request["params"] = params;
        
        QString json = QJsonDocument(request).toJson(QJsonDocument::Compact) + "\n";
        m_process->write(json.toUtf8());
        
        UQFF_LOG_DEBUG("PythonBridge", QString("Called: %1").arg(method));
    }
    
    // Compute UQFF equation
    void computeUQFF(int equationNum, const QJsonObject& systemParams, double t = 1.0) {
        QJsonObject params;
        params["equation"] = equationNum;
        params["system"] = systemParams;
        params["t"] = t;
        callAsync("compute_equation", params);
    }
    
    // Run validation tests
    void runTests(const QString& className = "") {
        QJsonObject params;
        params["class_name"] = className;
        callAsync("run_tests", params);
    }
    
    // Query astronomical system
    void querySystem(const QString& systemName) {
        QJsonObject params;
        params["query"] = systemName;
        callAsync("query_system", params);
    }
    
signals:
    void resultReady(const QString& method, const QJsonObject& result);
    void errorOccurred(const QString& error);
    void processFinished();
    
private slots:
    void onStdoutReady() {
        while (m_process->canReadLine()) {
            QString line = QString::fromUtf8(m_process->readLine()).trimmed();
            
            // Check for UQFF_LOG format from Python
            if (line.startsWith("UQFF_LOG:")) {
                QString jsonStr = line.mid(9);
                QJsonDocument doc = QJsonDocument::fromJson(jsonStr.toUtf8());
                if (doc.isObject()) {
                    QJsonObject logObj = doc.object();
                    LogLevel level = LogLevel::INFO;
                    QString levelStr = logObj["level"].toString();
                    if (levelStr == "DEBUG") level = LogLevel::DEBUG;
                    else if (levelStr == "PHYSICS") level = LogLevel::PHYSICS;
                    else if (levelStr == "VALIDATION") level = LogLevel::VALIDATION;
                    else if (levelStr == "WARNING") level = LogLevel::WARNING;
                    else if (levelStr == "ERROR") level = LogLevel::ERROR;
                    
                    UQFFEventBus::instance().log(
                        level,
                        logObj["source"].toString(),
                        logObj["message"].toString(),
                        logObj["data"].toObject()
                    );
                }
                continue;
            }
            
            // Try to parse as JSON-RPC response
            QJsonDocument doc = QJsonDocument::fromJson(line.toUtf8());
            if (doc.isObject()) {
                QJsonObject response = doc.object();
                QString method = response["method"].toString();
                QJsonObject result = response["result"].toObject();
                
                emit resultReady(method, result);
                
                // Also emit to event bus
                if (result.contains("equation")) {
                    UQFFEventBus::instance().emitComputationCompleted(
                        "Python", result["equation"].toString(), result);
                }
            }
        }
    }
    
    void onStderrReady() {
        QString error = QString::fromUtf8(m_process->readAllStandardError());
        if (!error.trimmed().isEmpty()) {
            UQFF_LOG_ERROR("Python", error.trimmed());
            emit errorOccurred(error);
        }
    }
    
    void onProcessFinished(int exitCode, QProcess::ExitStatus status) {
        m_running = false;
        if (status == QProcess::CrashExit) {
            UQFF_LOG_ERROR("PythonBridge", QString("Python process crashed (code %1)").arg(exitCode));
        } else {
            UQFF_LOG_INFO("PythonBridge", QString("Python process exited (code %1)").arg(exitCode));
        }
        emit processFinished();
    }
    
private:
    QString findPythonExecutable() {
        // Check for venv first
        QString venvPython = ".venv/Scripts/python.exe";
        if (QFile::exists(venvPython)) {
            return venvPython;
        }
        // Fall back to system Python
        return "python";
    }
    
    QProcess* m_process;
    bool m_running = false;
    int m_requestId = 1;
    QString m_scriptPath = "CondensedPhysics.py";
    
public:
    void setPythonScript(const QString& path) { m_scriptPath = path; }
    QString pythonScript() const { return m_scriptPath; }
};


// ═══════════════════════════════════════════════════════════════════════════════
// 3. COMPARISON DASHBOARD - Split-View Results Comparison (Tab 10)
// ═══════════════════════════════════════════════════════════════════════════════

class ComparisonDashboard : public QWidget {
    Q_OBJECT
    
public:
    ComparisonDashboard(QWidget* parent = nullptr) : QWidget(parent) {
        setupUI();
        connectEventBus();
        
        UQFF_LOG_INFO("ComparisonDashboard", "Dashboard initialized");
    }
    
private:
    void setupUI() {
        QVBoxLayout* mainLayout = new QVBoxLayout(this);
        
        // Header
        QLabel* header = new QLabel(
            "<h2 style='color: #9C27B0;'>📊 Comparison Dashboard</h2>"
            "<p style='color: #888;'>Real-time comparison: MAIN_1_CoAnQi.cpp vs CondensedPhysics.py</p>"
        );
        mainLayout->addWidget(header);
        
        // System selector
        QHBoxLayout* selectorLayout = new QHBoxLayout();
        selectorLayout->addWidget(new QLabel("System:"));
        systemSelector = new QComboBox();
        systemSelector->addItems({"Sun", "Sagittarius A*", "Betelgeuse", "NGC 3596", "Custom..."});
        systemSelector->setEditable(true);
        selectorLayout->addWidget(systemSelector);
        
        QPushButton* compareBtn = new QPushButton("🔄 Compare");
        connect(compareBtn, &QPushButton::clicked, this, &ComparisonDashboard::runComparison);
        selectorLayout->addWidget(compareBtn);
        
        selectorLayout->addStretch();
        mainLayout->addLayout(selectorLayout);
        
        // Split comparison table
        QSplitter* splitter = new QSplitter(Qt::Horizontal);
        
        // Left side: C++ results
        QGroupBox* cppGroup = new QGroupBox("🔧 MAIN_1_CoAnQi.cpp (C++)");
        QVBoxLayout* cppLayout = new QVBoxLayout(cppGroup);
        cppTable = createResultsTable();
        cppLayout->addWidget(cppTable);
        splitter->addWidget(cppGroup);
        
        // Right side: Python results
        QGroupBox* pyGroup = new QGroupBox("🐍 CondensedPhysics.py (Python)");
        QVBoxLayout* pyLayout = new QVBoxLayout(pyGroup);
        pythonTable = createResultsTable();
        pyLayout->addWidget(pythonTable);
        splitter->addWidget(pyGroup);
        
        mainLayout->addWidget(splitter, 1);
        
        // Difference summary
        QGroupBox* diffGroup = new QGroupBox("📈 Difference Analysis");
        QVBoxLayout* diffLayout = new QVBoxLayout(diffGroup);
        
        diffTable = new QTableWidget();
        diffTable->setColumnCount(5);
        diffTable->setHorizontalHeaderLabels({"Component", "C++ Value", "Python Value", "Δ (Abs)", "Δ (%)"});
        diffTable->horizontalHeader()->setStretchLastSection(true);
        diffTable->setStyleSheet("QTableWidget { font-family: 'Consolas', monospace; }");
        diffLayout->addWidget(diffTable);
        
        summaryLabel = new QLabel("Run a comparison to see results.");
        summaryLabel->setStyleSheet("color: #888; padding: 10px;");
        diffLayout->addWidget(summaryLabel);
        
        mainLayout->addWidget(diffGroup);
        
        setLayout(mainLayout);
    }
    
    QTableWidget* createResultsTable() {
        QTableWidget* table = new QTableWidget();
        table->setColumnCount(2);
        table->setHorizontalHeaderLabels({"Component", "Value (J/m³)"});
        table->horizontalHeader()->setStretchLastSection(true);
        table->setStyleSheet(
            "QTableWidget { font-family: 'Consolas', monospace; background: #1a1a1a; color: #e0e0e0; }"
            "QHeaderView::section { background: #333; color: #fff; padding: 5px; }"
        );
        
        // Pre-populate rows
        QStringList components = {"Ug1", "Ug2", "Ug3", "Ug4", "Ug_total", "U_b1", "Um", "F_U"};
        table->setRowCount(components.size());
        for (int i = 0; i < components.size(); i++) {
            table->setItem(i, 0, new QTableWidgetItem(components[i]));
            table->setItem(i, 1, new QTableWidgetItem("--"));
        }
        
        return table;
    }
    
    void connectEventBus() {
        connect(&UQFFEventBus::instance(), &UQFFEventBus::computationCompleted,
                this, &ComparisonDashboard::onComputationResult);
    }
    
private slots:
    void runComparison() {
        QString system = systemSelector->currentText();
        UQFF_LOG_INFO("ComparisonDashboard", QString("Running comparison for: %1").arg(system));
        
        m_cppResults.clear();
        m_pythonResults.clear();
        
        // Request computation from both sources via event bus
        QJsonObject params;
        params["system"] = system;
        params["equations"] = QJsonArray({1, 2, 3, 4, 5, 6, 7, 8});
        
        UQFFEventBus::instance().requestComputation("C++", "compute_all", params);
        UQFFEventBus::instance().requestComputation("Python", "compute_all", params);
    }
    
    void onComputationResult(const QString& source, const QString& equation, 
                              const QJsonObject& results) {
        // Store results based on source
        if (source.contains("C++") || source.contains("MAIN_1")) {
            m_cppResults[equation] = results;
            updateTable(cppTable, results);
        } else if (source.contains("Python") || source.contains("Condensed")) {
            m_pythonResults[equation] = results;
            updateTable(pythonTable, results);
        }
        
        // If we have both, compute differences
        if (!m_cppResults.isEmpty() && !m_pythonResults.isEmpty()) {
            computeDifferences();
        }
    }
    
    void updateTable(QTableWidget* table, const QJsonObject& results) {
        QStringList keys = {"Ug1", "Ug2", "Ug3", "Ug4", "Ug_total", "U_b1", "Um", "F_U"};
        for (int i = 0; i < keys.size() && i < table->rowCount(); i++) {
            if (results.contains(keys[i])) {
                double value = results[keys[i]].toDouble();
                table->item(i, 1)->setText(QString::number(value, 'e', 4));
            }
        }
    }
    
    void computeDifferences() {
        diffTable->setRowCount(0);
        
        QStringList keys = {"Ug1", "Ug2", "Ug3", "Ug4", "Ug_total", "U_b1", "Um", "F_U"};
        double maxError = 0;
        int passCount = 0;
        
        for (const QString& key : keys) {
            // Find values from both sources
            double cppVal = 0, pyVal = 0;
            for (auto& results : m_cppResults) {
                if (results.contains(key)) {
                    cppVal = results[key].toDouble();
                    break;
                }
            }
            for (auto& results : m_pythonResults) {
                if (results.contains(key)) {
                    pyVal = results[key].toDouble();
                    break;
                }
            }
            
            if (cppVal == 0 && pyVal == 0) continue;
            
            double absDiff = std::abs(cppVal - pyVal);
            double pctDiff = (cppVal != 0) ? (absDiff / std::abs(cppVal)) * 100 : 0;
            
            int row = diffTable->rowCount();
            diffTable->insertRow(row);
            
            diffTable->setItem(row, 0, new QTableWidgetItem(key));
            diffTable->setItem(row, 1, new QTableWidgetItem(QString::number(cppVal, 'e', 4)));
            diffTable->setItem(row, 2, new QTableWidgetItem(QString::number(pyVal, 'e', 4)));
            diffTable->setItem(row, 3, new QTableWidgetItem(QString::number(absDiff, 'e', 4)));
            diffTable->setItem(row, 4, new QTableWidgetItem(QString::number(pctDiff, 'f', 2) + "%"));
            
            // Color coding
            QString color = (pctDiff < 1) ? "#4CAF50" : (pctDiff < 10) ? "#FF9800" : "#F44336";
            for (int c = 0; c < 5; c++) {
                diffTable->item(row, c)->setBackground(QColor(color));
            }
            
            maxError = std::max(maxError, pctDiff);
            if (pctDiff < 1) passCount++;
        }
        
        // Summary
        QString status = (maxError < 1) ? "✅ EXCELLENT" : (maxError < 10) ? "⚠️ ACCEPTABLE" : "❌ DIVERGENT";
        summaryLabel->setText(QString("%1 | Max error: %2% | %3/%4 components within 1%")
                              .arg(status).arg(maxError, 0, 'f', 2)
                              .arg(passCount).arg(keys.size()));
        
        // Log validation result
        UQFFEventBus::instance().emitValidationResult(
            "C++ vs Python Comparison",
            maxError < 10,
            maxError,
            QString("%1/%2 within tolerance").arg(passCount).arg(keys.size())
        );
    }
    
private:
    QComboBox* systemSelector;
    QTableWidget* cppTable;
    QTableWidget* pythonTable;
    QTableWidget* diffTable;
    QLabel* summaryLabel;
    
    QHash<QString, QJsonObject> m_cppResults;
    QHash<QString, QJsonObject> m_pythonResults;
};


// ═══════════════════════════════════════════════════════════════════════════════
// 4. SESSION PERSISTENCE - Save/Restore Session State
// ═══════════════════════════════════════════════════════════════════════════════

class SessionPersistence : public QObject {
    Q_OBJECT
    
public:
    static SessionPersistence& instance() {
        static SessionPersistence sp;
        return sp;
    }
    
    // ═══════════════════════════════════════════════════════════════════════════
    // SAVE SESSION
    // ═══════════════════════════════════════════════════════════════════════════
    
    bool saveSession(const QString& filepath = "") {
        QString path = filepath;
        if (path.isEmpty()) {
            path = QString("sessions/uqff_session_%1.json")
                   .arg(QDateTime::currentDateTime().toString("yyyyMMdd_hhmmss"));
        }
        
        QJsonObject session;
        session["version"] = "1.0";
        session["timestamp"] = QDateTime::currentDateTime().toString(Qt::ISODate);
        session["application"] = "Source2 UQFF Calculator";
        
        // Save systems
        session["systems"] = m_systems;
        
        // Save parameters
        session["parameters"] = m_parameters;
        
        // Save simulation history (last 100)
        QJsonArray simHistory;
        int start = (m_simulationHistory.size() > 100) ? static_cast<int>(m_simulationHistory.size()) - 100 : 0;
        for (int i = start; i < static_cast<int>(m_simulationHistory.size()); i++) {
            simHistory.append(m_simulationHistory[i]);
        }
        session["simulation_history"] = simHistory;
        
        // Save log entries (last 500)
        QJsonArray logs;
        auto logHistory = UQFFEventBus::instance().getLogHistory();
        start = (logHistory.size() > 500) ? static_cast<int>(logHistory.size()) - 500 : 0;
        for (int i = start; i < static_cast<int>(logHistory.size()); i++) {
            logs.append(logHistory[i].toJson());
        }
        session["log_entries"] = logs;
        
        // Ensure directory exists
        QFileInfo fi(path);
        QDir().mkpath(fi.absolutePath());
        
        // Write file
        QFile file(path);
        if (file.open(QIODevice::WriteOnly)) {
            file.write(QJsonDocument(session).toJson(QJsonDocument::Indented));
            file.close();
            UQFF_LOG_INFO("Session", QString("Session saved: %1").arg(path));
            return true;
        }
        
        UQFF_LOG_ERROR("Session", QString("Failed to save session: %1").arg(path));
        return false;
    }
    
    // ═══════════════════════════════════════════════════════════════════════════
    // LOAD SESSION
    // ═══════════════════════════════════════════════════════════════════════════
    
    bool loadSession(const QString& filepath) {
        QFile file(filepath);
        if (!file.open(QIODevice::ReadOnly)) {
            UQFF_LOG_ERROR("Session", QString("Failed to open session: %1").arg(filepath));
            return false;
        }
        
        QJsonDocument doc = QJsonDocument::fromJson(file.readAll());
        file.close();
        
        if (!doc.isObject()) {
            UQFF_LOG_ERROR("Session", "Invalid session file format");
            return false;
        }
        
        QJsonObject session = doc.object();
        
        // Restore systems
        m_systems = session["systems"].toArray();
        
        // Restore parameters
        m_parameters = session["parameters"].toObject();
        
        // Restore simulation history
        m_simulationHistory.clear();
        for (const auto& item : session["simulation_history"].toArray()) {
            m_simulationHistory.append(item.toObject());
        }
        
        UQFF_LOG_INFO("Session", QString("Session loaded: %1 (%2 systems, %3 simulations)")
                      .arg(filepath)
                      .arg(m_systems.size())
                      .arg(m_simulationHistory.size()));
        
        emit sessionLoaded(session);
        return true;
    }
    
    // ═══════════════════════════════════════════════════════════════════════════
    // STATE MANAGEMENT
    // ═══════════════════════════════════════════════════════════════════════════
    
    void addSystem(const QJsonObject& system) {
        m_systems.append(system);
    }
    
    void setParameter(const QString& key, const QJsonValue& value) {
        m_parameters[key] = value;
    }
    
    QJsonValue getParameter(const QString& key) const {
        return m_parameters[key];
    }
    
    void recordSimulation(const QJsonObject& result) {
        QJsonObject entry;
        entry["timestamp"] = QDateTime::currentDateTime().toString(Qt::ISODate);
        entry["result"] = result;
        m_simulationHistory.append(entry);
    }
    
    QJsonArray getSystems() const { return m_systems; }
    QJsonObject getParameters() const { return m_parameters; }
    QList<QJsonObject> getSimulationHistory() const { return m_simulationHistory; }
    
signals:
    void sessionLoaded(const QJsonObject& session);
    void sessionSaved(const QString& filepath);
    
private:
    SessionPersistence() : QObject(nullptr) {}
    
    QJsonArray m_systems;
    QJsonObject m_parameters;
    QList<QJsonObject> m_simulationHistory;
};

// ═══════════════════════════════════════════════════════════════════════════════
// 5. UQFF JAVASCRIPT SERVER WIDGET - HTTP Client for uqff_server.js (Gap #5 Fix)
// ═══════════════════════════════════════════════════════════════════════════════

#include <QNetworkAccessManager>
#include <QNetworkRequest>
#include <QNetworkReply>
#include <QUrl>
#include <QUrlQuery>
#include <QProgressBar>

/**
 * @class UQFFJavaScriptWidget
 * @brief HTTP client widget for communicating with uqff_server.js (port 3141)
 * 
 * Bridges the Qt/C++ GUI to the 23,790-line index.js JavaScript engine
 * via the REST API provided by uqff_server.js.
 * 
 * Endpoints used:
 *   GET  /api/health     - Check if server is running
 *   GET  /api/systems    - List available systems
 *   POST /api/compute    - Compute UQFF for a system
 *   GET  /api/constants  - Get physics constants
 */
class UQFFJavaScriptWidget : public QWidget {
    Q_OBJECT

public:
    explicit UQFFJavaScriptWidget(QWidget* parent = nullptr)
        : QWidget(parent), serverAvailable(false) {
        
        setupUI();
        networkManager = new QNetworkAccessManager(this);
        
        connect(networkManager, &QNetworkAccessManager::finished,
                this, &UQFFJavaScriptWidget::handleNetworkReply);
        
        // Check server health on startup
        QTimer::singleShot(1000, this, &UQFFJavaScriptWidget::checkServerHealth);
    }

private:
    void setupUI() {
        QVBoxLayout* mainLayout = new QVBoxLayout(this);
        mainLayout->setSpacing(10);
        mainLayout->setContentsMargins(15, 15, 15, 15);
        
        // ═══════════════════════════════════════════════════════════════════════
        // HEADER
        // ═══════════════════════════════════════════════════════════════════════
        
        QLabel* header = new QLabel(this);
        header->setText(
            "<h2 style='color: #FF9800;'>🌐 UQFF JavaScript Engine</h2>"
            "<p style='color: #888;'>HTTP client for uqff_server.js (port 3141)</p>"
        );
        mainLayout->addWidget(header);
        
        // ═══════════════════════════════════════════════════════════════════════
        // SERVER STATUS
        // ═══════════════════════════════════════════════════════════════════════
        
        QHBoxLayout* statusLayout = new QHBoxLayout();
        
        statusLabel = new QLabel("⏳ Checking server...", this);
        statusLabel->setStyleSheet("font-size: 11pt; padding: 5px;");
        statusLayout->addWidget(statusLabel);
        
        statusLayout->addStretch();
        
        QPushButton* startServerBtn = new QPushButton("▶ Start Server", this);
        startServerBtn->setStyleSheet(
            "background-color: #4CAF50; "
            "color: white; "
            "padding: 8px 16px; "
            "font-weight: bold;"
        );
        connect(startServerBtn, &QPushButton::clicked, this, &UQFFJavaScriptWidget::startServer);
        statusLayout->addWidget(startServerBtn);
        
        QPushButton* refreshBtn = new QPushButton("🔄", this);
        refreshBtn->setToolTip("Check server status");
        connect(refreshBtn, &QPushButton::clicked, this, &UQFFJavaScriptWidget::checkServerHealth);
        statusLayout->addWidget(refreshBtn);
        
        mainLayout->addLayout(statusLayout);
        
        // ═══════════════════════════════════════════════════════════════════════
        // COMPUTATION INPUT
        // ═══════════════════════════════════════════════════════════════════════
        
        QGroupBox* computeGroup = new QGroupBox("UQFF Computation", this);
        QVBoxLayout* computeLayout = new QVBoxLayout(computeGroup);
        
        // System selector
        QHBoxLayout* systemLayout = new QHBoxLayout();
        systemLayout->addWidget(new QLabel("System:", this));
        
        systemCombo = new QComboBox(this);
        systemCombo->addItem("Sagittarius A* (SgrA*)");
        systemCombo->addItem("SGR1745 Magnetar");
        systemCombo->addItem("M87* Black Hole");
        systemCombo->addItem("Sun");
        systemCombo->addItem("NGC 3596");
        systemCombo->addItem("Custom");
        systemLayout->addWidget(systemCombo);
        
        computeLayout->addLayout(systemLayout);
        
        // Custom parameters
        QHBoxLayout* paramsLayout = new QHBoxLayout();
        
        paramsLayout->addWidget(new QLabel("M:", this));
        massInput = new QLineEdit("8.155e36", this);
        massInput->setMaximumWidth(120);
        paramsLayout->addWidget(massInput);
        
        paramsLayout->addWidget(new QLabel("r:", this));
        radiusInput = new QLineEdit("4.4e19", this);
        radiusInput->setMaximumWidth(120);
        paramsLayout->addWidget(radiusInput);
        
        paramsLayout->addWidget(new QLabel("B₀:", this));
        bFieldInput = new QLineEdit("1e-4", this);
        bFieldInput->setMaximumWidth(120);
        paramsLayout->addWidget(bFieldInput);
        
        paramsLayout->addStretch();
        computeLayout->addLayout(paramsLayout);
        
        // Compute button
        QPushButton* computeBtn = new QPushButton("🔬 Compute UQFF", this);
        computeBtn->setStyleSheet(
            "background-color: #2196F3; "
            "color: white; "
            "padding: 10px 20px; "
            "font-size: 12pt; "
            "font-weight: bold;"
        );
        connect(computeBtn, &QPushButton::clicked, this, &UQFFJavaScriptWidget::computeUQFF);
        computeLayout->addWidget(computeBtn);
        
        mainLayout->addWidget(computeGroup);
        
        // Progress bar
        progressBar = new QProgressBar(this);
        progressBar->setVisible(false);
        mainLayout->addWidget(progressBar);
        
        // ═══════════════════════════════════════════════════════════════════════
        // RESULTS DISPLAY
        // ═══════════════════════════════════════════════════════════════════════
        
        resultsDisplay = new QTextEdit(this);
        resultsDisplay->setReadOnly(true);
        resultsDisplay->setStyleSheet(
            "background-color: #1E1E1E; "
            "color: #D4D4D4; "
            "font-family: 'Consolas', 'Courier New', monospace; "
            "font-size: 10pt; "
            "padding: 10px;"
        );
        resultsDisplay->setPlainText(
            "═══════════════════════════════════════════════════════════════════════════════\n"
            " UQFF JavaScript Engine (index.js)\n"
            " Connected via: uqff_server.js REST API (port 3141)\n"
            "═══════════════════════════════════════════════════════════════════════════════\n\n"
            " This widget provides HTTP client access to the 23,790-line index.js\n"
            " JavaScript engine, which contains 106 astrophysical systems.\n\n"
            " To start the server:\n"
            "   1. Open a terminal\n"
            "   2. Run: node uqff_server.js\n"
            "   3. Server will listen on http://127.0.0.1:3141\n\n"
            " Or click 'Start Server' above to auto-launch.\n"
            "═══════════════════════════════════════════════════════════════════════════════\n"
        );
        mainLayout->addWidget(resultsDisplay);
        
        setLayout(mainLayout);
    }

private slots:
    void checkServerHealth() {
        QUrl url("http://127.0.0.1:3141/api/health");
        QNetworkRequest request(url);
        request.setAttribute(QNetworkRequest::User, QString("health"));
        
        statusLabel->setText("⏳ Checking server...");
        networkManager->get(request);
    }
    
    void startServer() {
        // Start uqff_server.js via QProcess
        serverProcess = new QProcess(this);
        serverProcess->setWorkingDirectory(QCoreApplication::applicationDirPath());
        
        // Try to find node.js executable
        QStringList nodeExes = {"node", "node.exe", "C:/Program Files/nodejs/node.exe"};
        QString nodeExe;
        for (const QString& exe : nodeExes) {
            if (QFile::exists(exe) || exe == "node") {
                nodeExe = exe;
                break;
            }
        }
        
        QString serverPath = QCoreApplication::applicationDirPath() + "/uqff_server.js";
        if (!QFile::exists(serverPath)) {
            serverPath = QCoreApplication::applicationDirPath() + "/../uqff_server.js";
        }
        
        if (QFile::exists(serverPath) || serverPath.contains("uqff_server.js")) {
            serverProcess->start(nodeExe, QStringList() << serverPath);
            
            if (serverProcess->waitForStarted(5000)) {
                resultsDisplay->append("\n✅ Server process started (PID: " + 
                                       QString::number(serverProcess->processId()) + ")\n");
                resultsDisplay->append("   Waiting for server to initialize...\n\n");
                
                // Check health after delay
                QTimer::singleShot(2000, this, &UQFFJavaScriptWidget::checkServerHealth);
            } else {
                resultsDisplay->append("\n❌ Failed to start server process\n");
                resultsDisplay->append("   Make sure Node.js is installed and in PATH\n\n");
            }
        } else {
            resultsDisplay->append("\n❌ uqff_server.js not found\n");
            resultsDisplay->append("   Expected: " + serverPath + "\n\n");
        }
    }
    
    void computeUQFF() {
        if (!serverAvailable) {
            resultsDisplay->append("\n❌ Server not available. Start the server first.\n");
            return;
        }
        
        progressBar->setVisible(true);
        progressBar->setRange(0, 0);  // Indeterminate
        
        // Build request JSON
        QJsonObject params;
        params["M"] = massInput->text().toDouble();
        params["r"] = radiusInput->text().toDouble();
        params["B0"] = bFieldInput->text().toDouble();
        params["system"] = systemCombo->currentText().split(" ").first();
        
        QJsonDocument doc(params);
        QByteArray jsonData = doc.toJson();
        
        QUrl url("http://127.0.0.1:3141/api/compute");
        QNetworkRequest request(url);
        request.setHeader(QNetworkRequest::ContentTypeHeader, "application/json");
        request.setAttribute(QNetworkRequest::User, QString("compute"));
        
        resultsDisplay->append("\n🔄 Computing UQFF for: " + systemCombo->currentText() + "\n");
        resultsDisplay->append("   M = " + massInput->text() + " kg\n");
        resultsDisplay->append("   r = " + radiusInput->text() + " m\n");
        resultsDisplay->append("   B₀ = " + bFieldInput->text() + " T\n\n");
        
        networkManager->post(request, jsonData);
    }
    
    void handleNetworkReply(QNetworkReply* reply) {
        progressBar->setVisible(false);
        
        QString requestType = reply->request().attribute(QNetworkRequest::User).toString();
        
        if (reply->error() != QNetworkReply::NoError) {
            if (requestType == "health") {
                serverAvailable = false;
                statusLabel->setText("❌ Server offline");
                statusLabel->setStyleSheet("color: #FF5252; font-size: 11pt; padding: 5px;");
            } else {
                resultsDisplay->append("❌ Error: " + reply->errorString() + "\n");
            }
            reply->deleteLater();
            return;
        }
        
        QByteArray responseData = reply->readAll();
        QJsonDocument doc = QJsonDocument::fromJson(responseData);
        QJsonObject response = doc.object();
        
        if (requestType == "health") {
            serverAvailable = true;
            statusLabel->setText("✅ Server online (port 3141)");
            statusLabel->setStyleSheet("color: #4CAF50; font-size: 11pt; padding: 5px;");
            
            resultsDisplay->append("\n✅ Connected to uqff_server.js\n");
            resultsDisplay->append("   Status: " + response["status"].toString() + "\n");
            resultsDisplay->append("   Timestamp: " + response["timestamp"].toString() + "\n\n");
            
        } else if (requestType == "compute") {
            resultsDisplay->append("═══════════════════════════════════════════════════════════════════════════════\n");
            resultsDisplay->append(" COMPUTATION RESULTS\n");
            resultsDisplay->append("═══════════════════════════════════════════════════════════════════════════════\n\n");
            
            QJsonObject result = response["result"].toObject();
            
            resultsDisplay->append(QString("   Ug1 = %1\n").arg(result["Ug1"].toDouble(), 0, 'e', 6));
            resultsDisplay->append(QString("   Ug2 = %1\n").arg(result["Ug2"].toDouble(), 0, 'e', 6));
            resultsDisplay->append(QString("   Ug3 = %1\n").arg(result["Ug3"].toDouble(), 0, 'e', 6));
            resultsDisplay->append(QString("   Ug4 = %1\n").arg(result["Ug4"].toDouble(), 0, 'e', 6));
            resultsDisplay->append(QString("\n   F_U_Bi_i = %1\n").arg(result["F_U_Bi_i"].toDouble(), 0, 'e', 6));
            resultsDisplay->append(QString("   compressed_g = %1\n").arg(result["compressed_g"].toDouble(), 0, 'e', 6));
            
            resultsDisplay->append("\n═══════════════════════════════════════════════════════════════════════════════\n\n");
        }
        
        reply->deleteLater();
    }

private:
    QNetworkAccessManager* networkManager;
    QProcess* serverProcess = nullptr;
    QLabel* statusLabel;
    QComboBox* systemCombo;
    QLineEdit* massInput;
    QLineEdit* radiusInput;
    QLineEdit* bFieldInput;
    QProgressBar* progressBar;
    QTextEdit* resultsDisplay;
    bool serverAvailable;
};


#endif // SOURCE2_WIDGETS_ENHANCED_H
