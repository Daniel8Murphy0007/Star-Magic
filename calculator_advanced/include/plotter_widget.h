/**
 * Plotter Widget (QCustomPlot Integration)
 * Thread b6d9bc22 Priority 2 - Iteration #31 Visualization
 * 
 * Interactive 2D/3D plotting for equation visualization
 */

#pragma once

#include <QWidget>
#include <qcustomplot.h>
#include <vector>
#include <string>
#include <map>

namespace CalculatorAdvanced {

/**
 * @brief Plot type enumeration
 */
enum class PlotType {
    Line2D,           // Standard x-y line plot
    Scatter2D,        // Scatter points
    Parametric2D,     // Parametric curve
    Contour,          // Contour plot for f(x,y)
    Heatmap,          // 2D heatmap
    Parametric3D,     // 3D parametric (requires colormap projection)
    Vector2D          // Vector field
};

/**
 * @brief Plot data structure
 */
struct PlotData {
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;  // Optional for 3D/heatmap
    std::string title;
    std::string x_label;
    std::string y_label;
    std::string legend;
    PlotType type;
};

/**
 * @brief Interactive plotting widget
 */
class PlotterWidget : public QWidget {
    Q_OBJECT
    
public:
    explicit PlotterWidget(QWidget* parent = nullptr);
    ~PlotterWidget() override;
    
    /**
     * @brief Plot 2D data
     */
    void plot2D(const PlotData& data);
    
    /**
     * @brief Plot parametric curve (x(t), y(t))
     */
    void plotParametric2D(
        const std::vector<double>& x,
        const std::vector<double>& y,
        const std::string& title = "Parametric Curve"
    );
    
    /**
     * @brief Plot 3D parametric curve (projected to 2D with color)
     */
    void plotParametric3D(
        const std::vector<double>& x,
        const std::vector<double>& y,
        const std::vector<double>& z,
        const std::string& title = "3D Parametric Curve"
    );
    
    /**
     * @brief Plot contour lines for f(x,y) = c
     */
    void plotContour(
        const std::vector<std::vector<double>>& grid,
        double x_min, double x_max,
        double y_min, double y_max,
        const std::string& title = "Contour Plot"
    );
    
    /**
     * @brief Plot heatmap
     */
    void plotHeatmap(
        const std::vector<std::vector<double>>& data,
        const std::string& title = "Heatmap",
        const std::string& colormap = "viridis"
    );
    
    /**
     * @brief Plot vector field (x, y, dx, dy)
     */
    void plotVectorField(
        const std::vector<double>& x,
        const std::vector<double>& y,
        const std::vector<double>& dx,
        const std::vector<double>& dy,
        const std::string& title = "Vector Field"
    );
    
    /**
     * @brief Add plot layer (multiple curves on same axes)
     */
    void addPlot(const PlotData& data);
    
    /**
     * @brief Clear all plots
     */
    void clearPlots();
    
    /**
     * @brief Set axis labels
     */
    void setAxisLabels(const std::string& x_label, const std::string& y_label);
    
    /**
     * @brief Set plot title
     */
    void setTitle(const std::string& title);
    
    /**
     * @brief Enable/disable legend
     */
    void setLegendVisible(bool visible);
    
    /**
     * @brief Enable/disable grid
     */
    void setGridVisible(bool visible);
    
    /**
     * @brief Set axis ranges
     */
    void setXRange(double min, double max);
    void setYRange(double min, double max);
    
    /**
     * @brief Auto-scale axes to fit data
     */
    void autoScale();
    
    /**
     * @brief Export to image
     */
    bool exportToImage(const std::string& filename, int width = 800, int height = 600);
    
    /**
     * @brief Export to PDF
     */
    bool exportToPDF(const std::string& filename, int width = 800, int height = 600);
    
    /**
     * @brief Enable interactive features (zoom, pan, data cursor)
     */
    void setInteractive(bool enabled);
    
public slots:
    void onDataCursorMoved(QMouseEvent* event);
    void onZoomChanged();
    
signals:
    void plotClicked(double x, double y);
    void rangeChanged(double x_min, double x_max, double y_min, double y_max);
    
private:
    QCustomPlot* plot_;
    std::vector<PlotData> plot_data_;
    
    // Helper methods
    QCPColorMap* createColorMap(const std::vector<std::vector<double>>& data);
    QCPCurve* createParametricCurve(const std::vector<double>& x, const std::vector<double>& y);
    void setupInteractions();
};

/**
 * @brief Example usage:
 * 
 * PlotterWidget* plotter = new PlotterWidget();
 * 
 * // 1. Simple 2D function plot: y = x^2
 * std::vector<double> x(1000);
 * std::vector<double> y(1000);
 * for (int i = 0; i < 1000; ++i) {
 *     x[i] = -10 + i * 0.02;
 *     y[i] = x[i] * x[i];
 * }
 * PlotData data{x, y, {}, "Parabola", "x", "y", "y=x²", PlotType::Line2D};
 * plotter->plot2D(data);
 * 
 * // 2. Parametric curve: circle
 * std::vector<double> x_circle(1000), y_circle(1000);
 * for (int i = 0; i < 1000; ++i) {
 *     double t = 2 * M_PI * i / 1000.0;
 *     x_circle[i] = cos(t);
 *     y_circle[i] = sin(t);
 * }
 * plotter->plotParametric2D(x_circle, y_circle, "Unit Circle");
 * 
 * // 3. 3D helix (projected with color)
 * std::vector<double> x_helix(1000), y_helix(1000), z_helix(1000);
 * for (int i = 0; i < 1000; ++i) {
 *     double t = 4 * M_PI * i / 1000.0;
 *     x_helix[i] = cos(t);
 *     y_helix[i] = sin(t);
 *     z_helix[i] = t;
 * }
 * plotter->plotParametric3D(x_helix, y_helix, z_helix, "Helix");
 * 
 * // 4. UQFF F_U_Bi_i vs radius
 * std::vector<double> r(100), F_U(100);
 * for (int i = 0; i < 100; ++i) {
 *     r[i] = 1e14 + i * 1e14;  // 1e14 to 1e16 m
 *     F_U[i] = compute_F_U_Bi_i(r[i]);  // UQFF calculation
 * }
 * PlotData uqff_data{r, F_U, {}, "UQFF Buoyancy Force", "Radius (m)", "Force (N)", "F_U_Bi_i", PlotType::Line2D};
 * plotter->plot2D(uqff_data);
 * 
 * // 5. Export to file
 * plotter->exportToImage("uqff_plot.png", 1920, 1080);
 * plotter->exportToPDF("uqff_plot.pdf", 800, 600);
 */

} // namespace CalculatorAdvanced
