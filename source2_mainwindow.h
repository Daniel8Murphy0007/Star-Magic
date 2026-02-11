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
#include "UQFFResultsWidget.h"

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

private:
    BrowserWindow **browserWindows;
    UQFFResultsWidget* uqffResultsWidget;
};

#endif // SOURCE2_MAINWINDOW_H
