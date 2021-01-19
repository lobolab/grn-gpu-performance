// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#include "concentplotwidget.h"

#include "Simulator/simulator.h"
#include "Simulator/simstate.h"
#include "../Private/colorsconfig.h"

#include <QLabel>
#include <QVBoxLayout>
#include <QSettings>
#include <QFileDialog>
#include <qwt_curve_fitter.h>
#include <qwt_plot_renderer.h>

namespace LoboLab {

ConcentPlotWidget::ConcentPlotWidget(QWidget *parent)
  : QWidget(parent), simulator_(NULL), axisHvalues_(NULL), axisVvalues_(NULL),
    axisXvalues_(NULL), indH_(0), indV_(0), xSize_(1000) {
  initPlots();
  
  saveImageAction_ = new QAction(tr("Save image..."), this);
  saveImageAction_->setIcon(QIcon(":/Images/famfamfam_silk_icons/picture_save.png"));
  connect(saveImageAction_, SIGNAL(triggered()), this, SLOT(saveImage()));

  setContextMenuPolicy(Qt::ActionsContextMenu);
  addAction(saveImageAction_); // add the action to the context menu
}

ConcentPlotWidget::~ConcentPlotWidget() {
  clear();
}

void ConcentPlotWidget::initPlots() {
  plotH_ = new QwtPlot();
  plotV_ = new QwtPlot();
  plotX_ = new QwtPlot();

  plotH_->setCanvasBackground(Qt::white);
  plotV_->setCanvasBackground(Qt::white);
  plotX_->setCanvasBackground(Qt::white);

  plotH_->setAxisScale(QwtPlot::yLeft, 0, 1);
  plotV_->setAxisScale(QwtPlot::yLeft, 0, 1);
  plotX_->setAxisScale(QwtPlot::yLeft, 0, 1);
  plotX_->setAxisScale(QwtPlot::xBottom, 0, xSize_);

  plotH_->enableAxis(QwtPlot::xBottom, false);
  plotV_->enableAxis(QwtPlot::xBottom, false);
  plotX_->enableAxis(QwtPlot::xBottom, false);
  
  plotH_->setAxisFont(QwtPlot::yLeft, QFont("Arial", 8));
  plotV_->setAxisFont(QwtPlot::yLeft, QFont("Arial", 8));
  plotX_->setAxisFont(QwtPlot::yLeft, QFont("Arial", 8));
  
  PlotMagnifierZeroLimited* magnifier = new PlotMagnifierZeroLimited((QwtPlotCanvas *) plotH_->canvas() );
  magnifier->setAxisEnabled(QwtPlot::xBottom, false);

  magnifier = new PlotMagnifierZeroLimited((QwtPlotCanvas *) plotV_->canvas() );
  magnifier->setAxisEnabled(QwtPlot::xBottom, false);

  magnifier = new PlotMagnifierZeroLimited((QwtPlotCanvas *) plotX_->canvas() );
  magnifier->setAxisEnabled(QwtPlot::xBottom, false);


  QLabel *hPlotLabel = new QLabel("Lateral product concentrations");
  hPlotLabel->setAlignment(Qt::AlignCenter);
  QLabel *vPlotLabel = new QLabel("Anterior-posterior product concentrations");
  vPlotLabel->setAlignment(Qt::AlignCenter);

  QLabel *xPlotLabel = new QLabel("Single point product concentrations through time");
  xPlotLabel->setAlignment(Qt::AlignCenter);

  QHBoxLayout *hTopLay = new QHBoxLayout();
  hTopLay->addWidget(hPlotLabel);
  hTopLay->addWidget(vPlotLabel);

  QHBoxLayout *hLay = new QHBoxLayout();
  hLay->addWidget(plotH_);
  hLay->addWidget(plotV_);

  QVBoxLayout *lay = new QVBoxLayout();
  lay->addLayout(hTopLay);
  lay->addLayout(hLay);
  lay->addWidget(xPlotLabel);
  lay->addWidget(plotX_);

  setLayout(lay);
}

void ConcentPlotWidget::clear() {
  plotH_->detachItems();
  curvesH_.clear();
  plotV_->detachItems();
  curvesV_.clear();
  plotX_->detachItems();
  curvesX_.clear();
  delete [] axisHvalues_;
  delete [] axisVvalues_;
  delete [] axisXvalues_;
  axisHvalues_ = NULL;
  axisVvalues_ = NULL;
  axisXvalues_ = NULL;

  for (int i=0; i<curvesHdata_.size(); ++i)
    delete [] curvesHdata_.at(i);

  curvesHdata_.clear();

  for (int i=0; i<curvesVdata_.size(); ++i)
    delete [] curvesVdata_.at(i);

  curvesVdata_.clear();

  for (int i=0; i<curvesXdata_.size(); ++i)
    delete [] curvesXdata_.at(i);

  curvesXdata_.clear();
}

void ConcentPlotWidget::setSimulator(const Simulator *sim) {
  clear();

  simulator_ = sim;

  int width = simulator_->domainSize().width();
  int height = simulator_->domainSize().height();

  indH_ = height/2;
  indV_ = width/2;

  axisHvalues_ = new double[width + 1]; // +1 for making nice first step
  for (int i=0; i<=width; ++i)
    axisHvalues_[i] = i;

  axisVvalues_ = new double[height + 1];
  for (int i=0; i<=height; ++i)
    axisVvalues_[i] = i;

  axisXvalues_ = new double[xSize_ + 1];
  for (int i=0; i<=xSize_; ++i)
    axisXvalues_[i] = i;


  for (int i=simulator_->nProducts()-1; i>=0; --i) {
    QColor color = ColorsConfig::colorsCells[i % ColorsConfig::nColors];
    QPen pen(color);
    if(i < 3)
    	pen.setWidthF(4);
    else
      pen.setWidthF(3);

    QwtPlotCurve *curveH = new QwtPlotCurve();
    curveH->setStyle(QwtPlotCurve::Lines);
    curveH->setPen(pen);
    curveH->setRenderHint(QwtPlotItem::RenderAntialiased);
    double *curveHdata = new double[width + 1];
    curveH->setRawSamples(axisHvalues_, curveHdata, width);
    curveH->attach(plotH_);
    curvesHdata_.append(curveHdata);
    curvesH_.append(curveH);

    QwtPlotCurve *curveV = new QwtPlotCurve();
    curveV->setStyle(QwtPlotCurve::Lines);

    curveV->setPen(pen);
    curveV->setRenderHint(QwtPlotItem::RenderAntialiased);
    double *curveVdata = new double[height + 1];
    curveV->setRawSamples(axisVvalues_, curveVdata, height);
    curveV->attach(plotV_);
    curvesVdata_.append(curveVdata);
    curvesV_.append(curveV);

    QwtPlotCurve *curveX = new QwtPlotCurve();
    curveX->setStyle(QwtPlotCurve::Lines);
    curveX->setPen(pen);
    curveX->setRenderHint(QwtPlotItem::RenderAntialiased);
    double *curveXdata = new double[xSize_ + 1];
    for (int j = 0; j <= xSize_; ++j)
      curveXdata[j] = 0;
    curveX->setRawSamples(axisXvalues_, curveXdata, xSize_);
    curveX->attach(plotX_);
    curvesXdata_.append(curveXdata);
    curvesX_.append(curveX);
  }
}

void ConcentPlotWidget::setIndH(int ind) {
  indH_ = ind;
  updateCurveHdata();
  plotH_->replot();
}

void ConcentPlotWidget::setIndV(int ind) {
  indV_ = ind;
  updateCurveVdata();
  plotV_->replot();
}

void ConcentPlotWidget::updatePlots() {
  updateCurveXdata();
  updateCurveHdata();
  updateCurveVdata();
  plotV_->replot();
  plotH_->replot();
  plotX_->replot();
}

void ConcentPlotWidget::updateCurveHdata() {
  const SimState &state = simulator_->simulatedState();
  int n = curvesHdata_.size();
  for (int i=0; i<n; ++i) {
    double *curveHdata = curvesHdata_.at(i);
    const Eigen::MatrixXd &conc = state.product(n-i-1);
    curveHdata[0] = conc(0, indH_);
    int width = state.size().width();
    for (int j=0; j<width; ++j) {
      curveHdata[j+1] = conc(j, indH_);
    }
  }
}

void ConcentPlotWidget::updateCurveVdata() {
  const SimState &state = simulator_->simulatedState();
  int n = curvesVdata_.size();
  for (int i=0; i<n; ++i) {
    double *curveVdata = curvesVdata_.at(i);
    const Eigen::MatrixXd &conc = state.product(n-i-1);
    curveVdata[0] = conc(indV_, 0);
    int height = state.size().height();
    for (int j=0; j<height; ++j) {
      curveVdata[j+1] = conc(indV_,j);
    }
  }
}

void ConcentPlotWidget::updateCurveXdata() {
  const SimState &state = simulator_->simulatedState();
  int n = curvesXdata_.size();
  for (int i=0; i<n; ++i) {
    double *curveXdata = curvesXdata_.at(i);
    const Eigen::MatrixXd &conc = state.product(n-i-1);
    //for(int j=0; j<xSize; ++j)
    //{
    //	curveXdata[j] = curveXdata[j+1];
    //}
    // A ring buffer would be even better, but it is not easy to implement
    // in qwt
    memmove(curveXdata, curveXdata + 1, xSize_*sizeof(double));

    curveXdata[xSize_] = conc(indV_, indH_);
  }
}

//void ConcentPlotWidget::axisPlotHChanged()
//{
//	QwtInterval axisInterval = plotH->axisInterval(QwtPlot::yLeft);
//	if(axisInterval.minValue() < 0)
//		plotH->setAxisScale(QwtPlot::yLeft, 0, axisInterval.maxValue() -
//			axisInterval.minValue());
//}
//
//void ConcentPlotWidget::axisPlotVChanged()
//{
//
//}

void ConcentPlotWidget::PlotMagnifierZeroLimited::rescale(double factor) {
  factor = qAbs( factor );
  if ( factor == 1.0 || factor == 0.0 )
    return;

  bool doReplot = false;
  QwtPlot* plt = plot();

  const bool autoReplot = plt->autoReplot();
  plt->setAutoReplot( false );

  for ( int axisId = 0; axisId < QwtPlot::axisCnt; axisId++ ) {
    const QwtScaleDiv &scaleDiv = plt->axisScaleDiv( axisId );
    if ( isAxisEnabled( axisId ) && !scaleDiv.isEmpty() ) {
      //const double center =
      //    scaleDiv->lowerBound() + scaleDiv->range() / 2;
      //const double width_2 = scaleDiv->range() / 2 * factor;

      //plt->setAxisScale( axisId, center - width_2, center + width_2 );
      plt->setAxisScale(axisId, scaleDiv.lowerBound() * factor,
                                scaleDiv.upperBound() * factor);
      doReplot = true;
    }
  }

  plt->setAutoReplot( autoReplot );

  if ( doReplot )
    plt->replot();
}

// private slot
void ConcentPlotWidget::saveImage() {
  QwtPlotCanvas *plotCanvas = dynamic_cast<QwtPlotCanvas*>(
                                        childAt(mapFromGlobal(QCursor::pos())));
  QwtPlot *plot;
  if (plotCanvas) {
    plot = plotCanvas->plot(); 

    QSettings settings;
    settings.beginGroup("ConcentPlotWidget");
    QString lastDir = settings.value("lastDir").toString();

    QString fileName = QFileDialog::getSaveFileName(this, "Save image",
                       lastDir, "Images (*.svg *.pdf *.ps *.png *.jpg *.bmp)");

    if (!fileName.isEmpty()) {
      QString newDir = QFileInfo(fileName).absolutePath();
      if (newDir != lastDir)
        settings.setValue("lastDir", newDir);

      QSizeF renderSize = plotCanvas->size();
      renderSize /= 5;
      //renderSize.scale(100, 100, Qt::KeepAspectRatio);
      QwtPlotRenderer renderer;
      renderer.renderDocument(plot, fileName, renderSize, 100);
    }
  }
}

}