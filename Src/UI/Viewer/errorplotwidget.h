// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#pragma once

#include <QWidget>
#include <QAction>
#include <qwt_plot_magnifier.h>
#include <qwt_plot.h>
#include <qwt_plot_curve.h>
#include <qwt_plot_canvas.h>
#include <qwt_scale_draw.h>

namespace LoboLab {

class Deme;
class DB;

class ErrorPlotWidget : public QWidget {
  Q_OBJECT

  public:
    explicit ErrorPlotWidget(QWidget *parent = NULL);
    ~ErrorPlotWidget();

    int plotTimeBest(int searchId, DB *db);
    int plotTimeAll(int searchId, DB *db);
    int plotTimeDeme(int demeId, DB *db);

    int plotGeneBest(int searchId, DB *db);
    int plotGeneAll(int searchId, DB *db);
    int plotGeneDeme(int demeId, DB *db);

    void setLinScale();
    void setLogScale();
    void clear();

    private slots:
    void saveImage();

  private:
    class PlotMagnifierBottomFixed : public QwtPlotMagnifier {
    public:
      explicit PlotMagnifierBottomFixed(QwtPlotCanvas * canvas)
        : QwtPlotMagnifier(canvas) {}

    protected:
      virtual void rescale(double factor);
  };

    class Power10ScaleDraw : public QwtScaleDraw {
    public:
      Power10ScaleDraw() {}
      virtual ~Power10ScaleDraw() {}

    protected:
      virtual QwtText label(double v) const;
  };

    void initPlot();
    int addPlotTime(const QString &cond, DB *db);
    int addPlotGene(const QString &cond, DB *db);
    void setSimpleErrorAxis();
    void setFullErrorAxis();
    void setSimpleCompAxis();
    void setFullCompAxis();
    void setGenerationAxis();
    void setTimeAxis();

    QwtPlot *plot_;
    QList<double *> valuesList;

    QAction *saveImageAction_;
};

}