// ═══════════════════════════════════════════════════════════════════════════════
// SOURCE2_EVENT_BUS.H - Inter-Component Communication System
// ═══════════════════════════════════════════════════════════════════════════════
// 
// Provides a singleton event bus for cross-widget communication in Source2.
// All widgets can emit and receive events without direct coupling.
//
// Usage:
//   // Emit event
//   UQFFEventBus::instance().emitPhysicsUpdate("UQFFSimulator", params);
//   
//   // Connect to event
//   connect(&UQFFEventBus::instance(), &UQFFEventBus::physicsUpdate,
//           this, &MyWidget::onPhysicsUpdate);
//
// © 2025-2026 Daniel T. Murphy - Star-Magic UQFF Framework
// ═══════════════════════════════════════════════════════════════════════════════

#ifndef SOURCE2_EVENT_BUS_H
#define SOURCE2_EVENT_BUS_H

#include <QObject>
#include <QString>
#include <QJsonObject>
#include <QJsonArray>
#include <QDateTime>
#include <QMutex>
#include <QMutexLocker>

// ═══════════════════════════════════════════════════════════════════════════════
// LOG LEVELS
// ═══════════════════════════════════════════════════════════════════════════════

// Windows defines ERROR as a macro in wingdi.h - must undefine it
#ifdef ERROR
#undef ERROR
#endif

enum class LogLevel {
    DEBUG = 0,      // Detailed debugging information
    INFO = 1,       // General information
    PHYSICS = 2,    // Physics computation results
    VALIDATION = 3, // Test/validation results
    WARNING = 4,    // Warnings
    ERROR = 5       // Errors
};

inline QString logLevelToString(LogLevel level) {
    switch (level) {
        case LogLevel::DEBUG: return "DEBUG";
        case LogLevel::INFO: return "INFO";
        case LogLevel::PHYSICS: return "PHYSICS";
        case LogLevel::VALIDATION: return "VALIDATION";
        case LogLevel::WARNING: return "WARNING";
        case LogLevel::ERROR: return "ERROR";
        default: return "UNKNOWN";
    }
}

inline QString logLevelToColor(LogLevel level) {
    switch (level) {
        case LogLevel::DEBUG: return "#888888";      // Gray
        case LogLevel::INFO: return "#2196F3";       // Blue
        case LogLevel::PHYSICS: return "#4CAF50";    // Green
        case LogLevel::VALIDATION: return "#9C27B0"; // Purple
        case LogLevel::WARNING: return "#FF9800";    // Orange
        case LogLevel::ERROR: return "#F44336";      // Red
        default: return "#FFFFFF";
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// LOG ENTRY STRUCTURE
// ═══════════════════════════════════════════════════════════════════════════════

struct LogEntry {
    QDateTime timestamp;
    LogLevel level;
    QString source;      // Widget/component name
    QString message;
    QJsonObject data;    // Optional structured data
    
    QJsonObject toJson() const {
        QJsonObject obj;
        obj["timestamp"] = timestamp.toString(Qt::ISODate);
        obj["level"] = logLevelToString(level);
        obj["source"] = source;
        obj["message"] = message;
        obj["data"] = data;
        return obj;
    }
    
    static LogEntry fromJson(const QJsonObject& obj) {
        LogEntry entry;
        entry.timestamp = QDateTime::fromString(obj["timestamp"].toString(), Qt::ISODate);
        entry.source = obj["source"].toString();
        entry.message = obj["message"].toString();
        entry.data = obj["data"].toObject();
        
        QString levelStr = obj["level"].toString();
        if (levelStr == "DEBUG") entry.level = LogLevel::DEBUG;
        else if (levelStr == "INFO") entry.level = LogLevel::INFO;
        else if (levelStr == "PHYSICS") entry.level = LogLevel::PHYSICS;
        else if (levelStr == "VALIDATION") entry.level = LogLevel::VALIDATION;
        else if (levelStr == "WARNING") entry.level = LogLevel::WARNING;
        else entry.level = LogLevel::ERROR;
        
        return entry;
    }
};

// ═══════════════════════════════════════════════════════════════════════════════
// UQFF EVENT BUS - Singleton Pattern
// ═══════════════════════════════════════════════════════════════════════════════

class UQFFEventBus : public QObject {
    Q_OBJECT
    
public:
    // Singleton accessor
    static UQFFEventBus& instance() {
        static UQFFEventBus bus;
        return bus;
    }
    
    // Prevent copying
    UQFFEventBus(const UQFFEventBus&) = delete;
    UQFFEventBus& operator=(const UQFFEventBus&) = delete;
    
    // ═══════════════════════════════════════════════════════════════════════════
    // LOG MANAGEMENT
    // ═══════════════════════════════════════════════════════════════════════════
    
    void log(LogLevel level, const QString& source, const QString& message, 
             const QJsonObject& data = QJsonObject()) {
        QMutexLocker locker(&m_mutex);
        
        LogEntry entry;
        entry.timestamp = QDateTime::currentDateTime();
        entry.level = level;
        entry.source = source;
        entry.message = message;
        entry.data = data;
        
        m_logHistory.append(entry);
        
        // Limit history size (keep last 10000 entries)
        if (m_logHistory.size() > 10000) {
            m_logHistory.removeFirst();
        }
        
        emit logMessage(entry);
    }
    
    QList<LogEntry> getLogHistory() const {
        QMutexLocker locker(&m_mutex);
        return m_logHistory;
    }
    
    void clearLogHistory() {
        QMutexLocker locker(&m_mutex);
        m_logHistory.clear();
    }
    
    // ═══════════════════════════════════════════════════════════════════════════
    // PHYSICS EVENT EMITTERS
    // ═══════════════════════════════════════════════════════════════════════════
    
    void emitSystemParametersChanged(const QString& systemName, const QJsonObject& params) {
        log(LogLevel::INFO, "EventBus", QString("System parameters changed: %1").arg(systemName), params);
        emit systemParametersChanged(systemName, params);
    }
    
    void emitComputationCompleted(const QString& source, const QString& equationName, 
                                   const QJsonObject& results) {
        QJsonObject data;
        data["equation"] = equationName;
        data["results"] = results;
        log(LogLevel::PHYSICS, source, QString("Computation completed: %1").arg(equationName), data);
        emit computationCompleted(source, equationName, results);
    }
    
    void emitSimulationFrame(int frame, double time, const QJsonObject& fieldState) {
        QJsonObject data;
        data["frame"] = frame;
        data["time"] = time;
        data["state"] = fieldState;
        // Don't log every frame (too noisy), just emit
        emit simulationFrameUpdated(frame, time, fieldState);
    }
    
    void emitValidationResult(const QString& testName, bool passed, double error, 
                               const QString& details = "") {
        QJsonObject data;
        data["passed"] = passed;
        data["error"] = error;
        data["details"] = details;
        log(LogLevel::VALIDATION, "Validation", 
            QString("%1: %2").arg(testName).arg(passed ? "PASS" : "FAIL"), data);
        emit validationResult(testName, passed, error, details);
    }
    
    // ═══════════════════════════════════════════════════════════════════════════
    // API/DATA EVENT EMITTERS
    // ═══════════════════════════════════════════════════════════════════════════
    
    void emitApiDataReceived(const QString& query, const QString& source, 
                              const QJsonObject& data) {
        QJsonObject logData;
        logData["query"] = query;
        logData["api_source"] = source;
        log(LogLevel::INFO, "API", QString("Data received for: %1").arg(query), logData);
        emit apiDataReceived(query, source, data);
    }
    
    void emitCsvGenerated(const QString& filepath, int rowCount) {
        QJsonObject data;
        data["filepath"] = filepath;
        data["rows"] = rowCount;
        log(LogLevel::INFO, "CSV", QString("Generated: %1 (%2 rows)").arg(filepath).arg(rowCount), data);
        emit csvGenerated(filepath, rowCount);
    }
    
    // ═══════════════════════════════════════════════════════════════════════════
    // CROSS-WIDGET COMMUNICATION
    // ═══════════════════════════════════════════════════════════════════════════
    
    void broadcast(const QString& channel, const QJsonObject& message) {
        log(LogLevel::DEBUG, "Broadcast", QString("Channel: %1").arg(channel), message);
        emit broadcastMessage(channel, message);
    }
    
    void requestComputation(const QString& target, const QString& method, 
                            const QJsonObject& params) {
        QJsonObject data;
        data["target"] = target;
        data["method"] = method;
        data["params"] = params;
        log(LogLevel::DEBUG, "Request", QString("Computation request to %1").arg(target), data);
        emit computationRequested(target, method, params);
    }
    
signals:
    // ═══════════════════════════════════════════════════════════════════════════
    // SIGNALS - Connect to these in your widgets
    // ═══════════════════════════════════════════════════════════════════════════
    
    // Logging
    void logMessage(const LogEntry& entry);
    
    // Physics events
    void systemParametersChanged(const QString& systemName, const QJsonObject& params);
    void computationCompleted(const QString& source, const QString& equationName, 
                               const QJsonObject& results);
    void simulationFrameUpdated(int frame, double time, const QJsonObject& fieldState);
    void validationResult(const QString& testName, bool passed, double error, 
                           const QString& details);
    
    // API/Data events
    void apiDataReceived(const QString& query, const QString& apiSource, 
                          const QJsonObject& data);
    void csvGenerated(const QString& filepath, int rowCount);
    
    // Cross-widget communication
    void broadcastMessage(const QString& channel, const QJsonObject& message);
    void computationRequested(const QString& target, const QString& method, 
                               const QJsonObject& params);
    
private:
    UQFFEventBus() : QObject(nullptr) {}
    
    mutable QMutex m_mutex;
    QList<LogEntry> m_logHistory;
};

// ═══════════════════════════════════════════════════════════════════════════════
// CONVENIENCE MACROS
// ═══════════════════════════════════════════════════════════════════════════════

#define UQFF_LOG_DEBUG(source, msg) \
    UQFFEventBus::instance().log(LogLevel::DEBUG, source, msg)

#define UQFF_LOG_INFO(source, msg) \
    UQFFEventBus::instance().log(LogLevel::INFO, source, msg)

#define UQFF_LOG_PHYSICS(source, msg, data) \
    UQFFEventBus::instance().log(LogLevel::PHYSICS, source, msg, data)

#define UQFF_LOG_VALIDATION(source, msg, data) \
    UQFFEventBus::instance().log(LogLevel::VALIDATION, source, msg, data)

#define UQFF_LOG_WARNING(source, msg) \
    UQFFEventBus::instance().log(LogLevel::WARNING, source, msg)

#define UQFF_LOG_ERROR(source, msg) \
    UQFFEventBus::instance().log(LogLevel::ERROR, source, msg)

#endif // SOURCE2_EVENT_BUS_H
