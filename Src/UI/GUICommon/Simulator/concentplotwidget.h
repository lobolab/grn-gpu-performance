// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#pragma once

#include <QWidget>
#include <QAction>
#include <qwt_plot_magnifier.h>
#include <qwt_plot.h>
#include <qwt_plot_curve.h>
#include <qwt_plot_canvas.h>

namespace LoboLab {

class Simulator;

class ConcentPlotWidget : public QWidget {
  Q_OBJECT

 public:
  ConcentPlotWidget(QWidget *parent = NULL);
  ~ConcentPlotWidget();

  void setSimulator(const Simulator *sim);
  void setIndH(int ind);
  void setIndV(int ind);
  inline int indH() {return indH_;};
  inline int indV() {return indV_;};

  void clear();

  // Bug workaround: in Qwt 6.1.0, the graphs get giant when dragging the window!
  virtual QSize	sizeHint () const {return QSize(10,10); }
  virtual QSize	minimumSizeHint () const {return QSize(10,10); }

 public slots:
  void updatePlots();

 private slots:
  void saveImage();

 private:

  class PlotMagnifierZeroLimited : public QwtPlotMagnifier {
   public:
    explicit PlotMagnifierZeroLimited(QwtPlotCanvas * canvas)
      : QwtPlotMagnifier(canvas) {}

   protected:
    virtual void rescale(double factor);
  };

  void initPlots();
  double *createCurveData(int size);
  void updateCurveHdata();
  void updateCurveVdata();
  void updateCurveXdata();

  const Simulator *simulator_;

  QwtPlot *plotH_;
  QwtPlot *plotV_;
  QwtPlot *plotX_;
  QList<QwtPlotCurve*> curvesH_;
  QList<QwtPlotCurve*> curvesV_;
  QList<QwtPlotCurve*> curvesX_;
  QList<double*> curvesHdata_;
  QList<double*> curvesVdata_;
  QList<double*> curvesXdata_;
  double *axisHvalues_;
  double *axisVvalues_;
  double *axisXvalues_;

  int indH_;
  int indV_;
  int xSize_;
  
  QAction *saveImageAction_;
};

} // namespace LoboLab
