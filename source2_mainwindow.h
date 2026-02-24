#ifndef SOURCE2_MAINWINDOW_H
#define SOURCE2_MAINWINDOW_H

#include <QMainWindow>
#include <QApplication>
#include <QTabWidget>
#include <QDockWidget>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QPushButton>
#include <QLineEdit>
#include <QTextEdit>
#include <QLabel>
#include <QIcon>
#include <QUrl>
#include <QWidget>
#include <QProcess>
#include <QSystemTrayIcon>  // System tray icon for notification area
#include <QMenu>            // Context menu for tray icon
#include "UQFFResultsWidget.h"

// UQFF Integration Components (Phase 2 - Feb 2026)
#include "csv_body_reader.h"  // For UQFF::CelestialBodyCSV

#ifdef _WIN32
#include <windows.h>
#include <shellapi.h>
#endif

// Forward declarations
class BrowserWindow;

// MainWindow class declaration for Qt MOC processing
// Implementation is in source2.cpp
class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow();
    ~MainWindow();

private slots:
    // UQFF Physics Integration
    void computeUQFF(const QString& systemName);
    void parseAndDisplayUQFFResults(const QString& jsonStr);
    
    // CSV Body Loading (Phase 2 integration - Feb 2026)
    void loadBodiesFromCSV(const QString& csvPath = QString());
    void onBodiesLoaded(const std::vector<UQFF::CelestialBodyCSV>& bodies);

private:
    BrowserWindow **browserWindows;
    UQFFResultsWidget* uqffResultsWidget;
    std::vector<UQFF::CelestialBodyCSV> loadedBodies;  // Cached bodies from CSV
    QSystemTrayIcon* trayIcon;  // System tray icon
    QMenu* trayMenu;            // Tray icon context menu
    
    void setupSystemTrayIcon();  // Initialize tray icon
};

#endif // SOURCE2_MAINWINDOW_H