#include "../include/plotter_widget.h"
#include <QVBoxLayout>
#include <QLabel>
#include <QPushButton>

PlotterWidget::PlotterWidget(QWidget* parent)
    : QWidget(parent)
{
    setupUI();
}

PlotterWidget::~PlotterWidget() {}

void PlotterWidget::setupUI() {
    QVBoxLayout* layout = new QVBoxLayout(this);
    
    plotWidget_ = new QCustomPlot(this);
    plotWidget_->setMinimumHeight(400);
    
    layout->addWidget(plotWidget_);
    
    // Control panel
    QHBoxLayout* controlLayout = new QHBoxLayout();
    
    QPushButton* clearBtn = new QPushButton("Clear", this);
    connect(clearBtn, &QPushButton::clicked, this, &PlotterWidget::clearPlot);
    
    controlLayout->addWidget(clearBtn);
    controlLayout->addStretch();
    
    layout->addLayout(controlLayout);
}

void PlotterWidget::plot2D(const PlotData& data) {
    if (data.xValues.empty() || data.yValues.empty()) {
        return;
    }
    
    QVector<double> x = QVector<double>::fromStdVector(data.xValues);
    QVector<double> y = QVector<double>::fromStdVector(data.yValues);
    
    int graphIndex = plotWidget_->graphCount();
    plotWidget_->addGraph();
    plotWidget_->graph(graphIndex)->setData(x, y);
    plotWidget_->graph(graphIndex)->setName(QString::fromStdString(data.title));
    
    // Auto-range
    plotWidget_->rescaleAxes();
    plotWidget_->replot();
}

void PlotterWidget::plotScatter(const PlotData& data) {
    if (data.xValues.empty() || data.yValues.empty()) {
        return;
    }
    
    QVector<double> x = QVector<double>::fromStdVector(data.xValues);
    QVector<double> y = QVector<double>::fromStdVector(data.yValues);
    
    int graphIndex = plotWidget_->graphCount();
    plotWidget_->addGraph();
    plotWidget_->graph(graphIndex)->setData(x, y);
    plotWidget_->graph(graphIndex)->setLineStyle(QCPGraph::lsNone);
    plotWidget_->graph(graphIndex)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, 5));
    
    plotWidget_->rescaleAxes();
    plotWidget_->replot();
}

void PlotterWidget::plotParametric(const PlotData& data) {
    plot2D(data); // Same as 2D for now
}

void PlotterWidget::plotPolar(const PlotData& data) {
    // TODO: Implement polar plotting
    plot2D(data);
}

void PlotterWidget::plotContour(const PlotData& data) {
    // TODO: Implement contour plotting
}

void PlotterWidget::plot3DSurface(const PlotData& data) {
    // TODO: Implement 3D surface plotting
}

void PlotterWidget::plotVectorField(const PlotData& data) {
    // TODO: Implement vector field plotting
}

void PlotterWidget::clearPlot() {
    plotWidget_->clearGraphs();
    plotWidget_->replot();
}

void PlotterWidget::exportPlot(const std::string& filename) {
    QString qfilename = QString::fromStdString(filename);
    
    if (qfilename.endsWith(".png")) {
        plotWidget_->savePng(qfilename);
    } else if (qfilename.endsWith(".pdf")) {
        plotWidget_->savePdf(qfilename);
    } else if (qfilename.endsWith(".jpg")) {
        plotWidget_->saveJpg(qfilename);
    }
}
