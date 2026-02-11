#ifndef UQFF_RESULTS_WIDGET_H
#define UQFF_RESULTS_WIDGET_H

#include <QWidget>
#include <QPushButton>
#include <QTextEdit>
#include <QVBoxLayout>
#include <QMessageBox>
#include <QFile>
#include <QTextStream>
#include <QDateTime>
#include <QDialog>
#include <QLabel>
#include <nlohmann/json.hpp>
#include <cmath>

// VTK includes for Phase 2 visualization
#ifndef NO_VTK
#include <vtkSmartPointer.h>
#include <vtkChartXY.h>
#include <vtkPlot.h>
#include <vtkTable.h>
#include <vtkDoubleArray.h>
#include <vtkContextView.h>
#include <vtkContextScene.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkAxis.h>
#include <vtkNew.h>
#endif

using json = nlohmann::json;

/**
 * @brief Widget for displaying UQFF physics computation results
 * 
 * This widget displays results from MAIN_1_CoAnQi.cpp calculations including:
 * - F_U_Bi_i (unified field force)
 * - g_compressed (compressed gravity)
 * - Ug1-4 components (magnetic dipole, charge-reactivity, string rotation, vacuum concentration)
 * - Ubi (buoyancy force)
 * 
 * Provides export to CSV and hooks for VTK visualization (Phase 2).
 */
class UQFFResultsWidget : public QWidget {
    Q_OBJECT

public:
    explicit UQFFResultsWidget(QWidget* parent = nullptr)
        : QWidget(parent) {
        QVBoxLayout* layout = new QVBoxLayout(this);

        rawDataDisplay = new QTextEdit(this);
        rawDataDisplay->setReadOnly(true);
        rawDataDisplay->setStyleSheet("background-color: #1e1e1e; color: #d4d4d4; font-family: 'Consolas', monospace;");
        layout->addWidget(rawDataDisplay);

        exportBtn = new QPushButton("📊 Export to CSV", this);
        connect(exportBtn, &QPushButton::clicked, this, &UQFFResultsWidget::exportResults);
        layout->addWidget(exportBtn);

        visualizeBtn = new QPushButton("📈 Visualize (VTK)", this);
        visualizeBtn->setEnabled(false);  // Phase 2 feature
        connect(visualizeBtn, &QPushButton::clicked, this, &UQFFResultsWidget::visualizeResults);
        layout->addWidget(visualizeBtn);

        setLayout(layout);
    }

    /**
     * @brief Update widget with new UQFF computation results
     * @param data JSON object containing physics results from CoAnQi_Wrapper.py
     */
    void setResults(const json& data) {
        currentData = data;

        QString output = "=== UQFF PHYSICS RESULTS ===\n\n";

        if (data.contains("result") && data["result"].contains("F_U_Bi_i")) {
            auto& result = data["result"];
            output += QString("F_U_Bi_i: %1 N\n").arg(QString::number(result["F_U_Bi_i"].get<double>(), 'e', 6));
            output += QString("g_compressed: %1 m/s²\n\n").arg(QString::number(result["g_compressed"].get<double>(), 'e', 6));

            output += "=== UQFF Components ===\n";
            output += QString("Ug1 (Magnetic Dipole): %1 m/s²\n").arg(QString::number(result["Ug1"].get<double>(), 'e', 6));
            output += QString("Ug2 (Charge-Reactivity): %1 m/s²\n").arg(QString::number(result["Ug2"].get<double>(), 'e', 6));
            output += QString("Ug3 (String Rotation): %1 m/s²\n").arg(QString::number(result["Ug3"].get<double>(), 'e', 6));
            output += QString("Ug4 (Vacuum Concentration): %1 m/s²\n").arg(QString::number(result["Ug4"].get<double>(), 'e', 6));
            output += QString("Ubi (Buoyancy): %1 m/s²\n\n").arg(QString::number(result["Ubi"].get<double>(), 'e', 6));

            if (result.contains("system_params")) {
                auto& params = result["system_params"];
                output += "=== System Parameters ===\n";
                output += QString("Mass (M): %1 kg\n").arg(QString::number(params["M"].get<double>(), 'e', 6));
                output += QString("Radius (r): %1 m\n").arg(QString::number(params["r"].get<double>(), 'e', 6));
                if (params.contains("L_X"))
                    output += QString("Luminosity (L_X): %1 W\n").arg(QString::number(params["L_X"].get<double>(), 'e', 6));
                if (params.contains("B0"))
                    output += QString("Magnetic Field (B0): %1 T\n").arg(QString::number(params["B0"].get<double>(), 'e', 6));
                if (params.contains("omega0"))
                    output += QString("Angular Velocity (ω₀): %1 rad/s\n").arg(QString::number(params["omega0"].get<double>(), 'e', 6));
                if (params.contains("v_cm"))
                    output += QString("Center of Mass Velocity (v): %1 m/s\n").arg(QString::number(params["v_cm"].get<double>(), 'e', 6));
                if (params.contains("T"))
                    output += QString("Temperature (T): %1 K\n").arg(QString::number(params["T"].get<double>(), 'e', 6));
            }
        } else {
            output += "No results available\n";
        }

        rawDataDisplay->setText(output);
        
        // Enable visualization button if data is valid (Phase 2)
        if (!data.empty() && data.contains("result")) {
            visualizeBtn->setEnabled(true);
        }
    }

private slots:
    /**
     * @brief Export current results to CSV with timestamp
     */
    void exportResults() {
        if (currentData.empty()) {
            QMessageBox::warning(this, "No Data", "No results to export");
            return;
        }

        QString timestamp = QDateTime::currentDateTime().toString("yyyyMMdd_HHmmss");
        QString filename = QString("uqff_results_%1.csv").arg(timestamp);

        QFile file(filename);
        if (file.open(QIODevice::WriteOnly | QIODevice::Text)) {
            QTextStream out(&file);
            out << "Parameter,Value,Unit\n";

            if (currentData.contains("result")) {
                auto& result = currentData["result"];
                out << "F_U_Bi_i," << result["F_U_Bi_i"].get<double>() << ",N\n";
                out << "g_compressed," << result["g_compressed"].get<double>() << ",m/s²\n";
                out << "Ug1," << result["Ug1"].get<double>() << ",m/s²\n";
                out << "Ug2," << result["Ug2"].get<double>() << ",m/s²\n";
                out << "Ug3," << result["Ug3"].get<double>() << ",m/s²\n";
                out << "Ug4," << result["Ug4"].get<double>() << ",m/s²\n";
                out << "Ubi," << result["Ubi"].get<double>() << ",m/s²\n";
            }

            file.close();
            QMessageBox::information(this, "Export Complete", QString("Results saved to %1").arg(filename));
        } else {
            QMessageBox::critical(this, "Export Failed", "Could not open file for writing");
        }
    }

    /**
     * @brief Phase 2: VTK visualization of UQFF force components
     * 
     * Creates interactive 3D visualization showing:
     * - Ug1 (Magnetic Dipole) - Red
     * - Ug2 (Charge-Reactivity) - Blue  
     * - Ug3 (String Rotation) - Green
     * - Ug4 (Vacuum Concentration) - Purple
     * - Ubi (Buoyancy) - Orange
     * - g_compressed (Total Gravity) - Black (thick line)
     */
    void visualizeResults() {
#ifdef NO_VTK
        QMessageBox::warning(this, "VTK Not Available", 
            "VTK visualization requires VTK library.\n"
            "Rebuild with VTK support to enable this feature.");
        return;
#else
        if (currentData.empty() || !currentData.contains("result")) {
            QMessageBox::warning(this, "No Data", "No results to visualize");
            return;
        }

        // Show info message before opening native VTK window
        QMessageBox::information(this, "Opening VTK Window",
            "Interactive VTK visualization will open in a new window.\n\n"
            "🔍 Controls:\n"
            "  • Left-click and drag to pan\n"
            "  • Mouse wheel to zoom\n"
            "  • Close window when done");

        // Create VTK Chart
        vtkNew<vtkChartXY> chart;
        chart->SetShowLegend(true);
        chart->SetTitle("UQFF Force Components vs Radius");
        chart->GetAxis(vtkAxis::BOTTOM)->SetTitle("Radius (normalized, log scale)");
        chart->GetAxis(vtkAxis::LEFT)->SetTitle("Force Magnitude (m/s², log scale)");

        // Create data table
        vtkNew<vtkTable> table;

        // Generate radius points (logarithmic scale from 0.1 to 100 system radii)
        vtkNew<vtkDoubleArray> radiusArray;
        radiusArray->SetName("Radius");
        
        int numPoints = 100;
        double r_system = 1.0;  // Normalized to system radius
        if (currentData["result"].contains("system_params") && 
            currentData["result"]["system_params"].contains("r")) {
            r_system = currentData["result"]["system_params"]["r"].get<double>();
        }

        // Create arrays for each force component
        vtkNew<vtkDoubleArray> ug1Array, ug2Array, ug3Array, ug4Array, ubiArray, gArray;
        ug1Array->SetName("Ug1 (Magnetic)");
        ug2Array->SetName("Ug2 (Charge)");
        ug3Array->SetName("Ug3 (String)");
        ug4Array->SetName("Ug4 (Vacuum)");
        ubiArray->SetName("Ubi (Buoyancy)");
        gArray->SetName("g_compressed (Total)");

        // Get current values
        auto& result = currentData["result"];
        double ug1_val = std::abs(result["Ug1"].get<double>());
        double ug2_val = std::abs(result["Ug2"].get<double>());
        double ug3_val = std::abs(result["Ug3"].get<double>());
        double ug4_val = std::abs(result["Ug4"].get<double>());
        double ubi_val = std::abs(result["Ubi"].get<double>());
        double g_val = std::abs(result["g_compressed"].get<double>());

        // Simulate force scaling with radius (inverse square + corrections)
        for (int i = 0; i < numPoints; ++i) {
            double r = std::pow(10.0, -1.0 + 3.0 * i / (numPoints - 1.0));  // 0.1 to 100 (normalized)
            double scale = 1.0 / (r * r);  // Inverse square law

            radiusArray->InsertNextValue(r);
            ug1Array->InsertNextValue(ug1_val * scale);
            ug2Array->InsertNextValue(ug2_val * scale);
            ug3Array->InsertNextValue(ug3_val * scale * (1.0 + 0.1 * std::sin(10.0 * r)));  // String oscillation
            ug4Array->InsertNextValue(ug4_val * scale);
            ubiArray->InsertNextValue(ubi_val * scale);
            gArray->InsertNextValue(g_val * scale);
        }

        // Add arrays to table
        table->AddColumn(radiusArray);
        table->AddColumn(ug1Array);
        table->AddColumn(ug2Array);
        table->AddColumn(ug3Array);
        table->AddColumn(ug4Array);
        table->AddColumn(ubiArray);
        table->AddColumn(gArray);

        // Add plots to chart with color coding
        vtkPlot* plotUg1 = chart->AddPlot(vtkChart::LINE);
        plotUg1->SetInputData(table, 0, 1);
        plotUg1->SetColor(255, 0, 0);  // Red - Magnetic
        plotUg1->SetWidth(2.0);

        vtkPlot* plotUg2 = chart->AddPlot(vtkChart::LINE);
        plotUg2->SetInputData(table, 0, 2);
        plotUg2->SetColor(0, 0, 255);  // Blue - Charge
        plotUg2->SetWidth(2.0);

        vtkPlot* plotUg3 = chart->AddPlot(vtkChart::LINE);
        plotUg3->SetInputData(table, 0, 3);
        plotUg3->SetColor(0, 200, 0);  // Green - String
        plotUg3->SetWidth(2.0);

        vtkPlot* plotUg4 = chart->AddPlot(vtkChart::LINE);
        plotUg4->SetInputData(table, 0, 4);
        plotUg4->SetColor(128, 0, 128);  // Purple - Vacuum
        plotUg4->SetWidth(2.0);

        vtkPlot* plotUbi = chart->AddPlot(vtkChart::LINE);
        plotUbi->SetInputData(table, 0, 5);
        plotUbi->SetColor(255, 128, 0);  // Orange - Buoyancy
        plotUbi->SetWidth(2.0);

        vtkPlot* plotG = chart->AddPlot(vtkChart::LINE);
        plotG->SetInputData(table, 0, 6);
        plotG->SetColor(0, 0, 0);  // Black - Total gravity
        plotG->SetWidth(3.5);  // Thicker line

        // Set logarithmic axes
        chart->GetAxis(vtkAxis::BOTTOM)->SetLogScale(true);
        chart->GetAxis(vtkAxis::LEFT)->SetLogScale(true);

        // Create context view for 2D charts
        vtkNew<vtkContextView> view;
        view->GetScene()->AddItem(chart);
        view->GetRenderWindow()->SetSize(1200, 800);
        view->GetRenderWindow()->SetWindowName("UQFF Force Field Visualization - Phase 2");

        // Start interaction (blocking until window closed)
        view->GetInteractor()->Initialize();
        view->GetInteractor()->Start();
#endif
    }

private:
    QTextEdit* rawDataDisplay;
    QPushButton* exportBtn;
    QPushButton* visualizeBtn;
    json currentData;
};

#endif // UQFF_RESULTS_WIDGET_H
