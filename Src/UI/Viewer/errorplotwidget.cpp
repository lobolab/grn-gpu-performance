// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#include "errorplotwidget.h"

#include "Search/deme.h"
#include "Search/generation.h"
#include "Search/individual.h"
#include "DB/db.h"

#include <QLabel>
#include <QVBoxLayout>
#include <QSettings>
#include <QFileDialog>
//#include <qwt_plot_marker.h>
#include <qwt_scale_engine.h>
#include <qwt_scale_engine.h>
#include <qwt_plot_renderer.h>


namespace LoboLab {

ErrorPlotWidget::ErrorPlotWidget(QWidget *parent)
  : QWidget(parent) {
  initPlot();

  saveImageAction_ = new QAction(tr("Save image..."), this);
  saveImageAction_->setIcon(QIcon(":/Images/famfamfam_silk_icons/picture_save.png"));
  connect(saveImageAction_, SIGNAL(triggered()), this, SLOT(saveImage()));

  setContextMenuPolicy(Qt::ActionsContextMenu);
  addAction(saveImageAction_); // add the action to the context menu
}

ErrorPlotWidget::~ErrorPlotWidget() {
  clear();
}

void ErrorPlotWidget::initPlot() {
  plot_ = new QwtPlot();
  
  plot_->setAxisScaleEngine(QwtPlot::xBottom, new QwtLogScaleEngine());
  plot_->setAxisScaleDraw(QwtPlot::xBottom, new Power10ScaleDraw());

  plot_->enableAxis(QwtPlot::yRight, true); // Complexity axis

  plot_->setAxisFont(QwtPlot::xBottom, QFont("Arial", 8));
  plot_->setAxisFont(QwtPlot::yLeft, QFont("Arial", 8));
  plot_->setAxisFont(QwtPlot::yRight, QFont("Arial", 8));

  plot_->setCanvasBackground(Qt::white);

  PlotMagnifierBottomFixed* magnifier =
    new PlotMagnifierBottomFixed((QwtPlotCanvas *)plot_->canvas());
  magnifier->setAxisEnabled(QwtPlot::xBottom, false);
  magnifier->setAxisEnabled(QwtPlot::yRight, false);

  QVBoxLayout *lay = new QVBoxLayout();
  lay->addWidget(plot_);
  lay->setMargin(0);
  setLayout(lay);
}

void ErrorPlotWidget::setSimpleErrorAxis() {
  QwtText text("<font color='blue'>Error</font>"); // (<font color='blue'>min</font>)");
  text.setFont(QFont("Arial", 10));
  plot_->setAxisTitle(QwtPlot::yLeft, text);
}

void ErrorPlotWidget::setFullErrorAxis() {
  QwtText text("Error (<font color='blue'>min</font>, "
    "<font color='green'>mean</font>, "
    "<font color='red'>max</font>)");
  text.setFont(QFont("Arial", 10));
  plot_->setAxisTitle(QwtPlot::yLeft, text);
}

void ErrorPlotWidget::setSimpleCompAxis() {
  QwtText text("Complexity"); // (<font color='#808080'>best</font>)");
  text.setFont(QFont("Arial", 10));
  plot_->setAxisTitle(QwtPlot::yRight, text);
}

void ErrorPlotWidget::setFullCompAxis() {
  QwtText text("Complexity (<font color='#808080'>best</font>, "
    "<font color='#c0c0c0'>max</font>, "
    "<font color='#c0c0c0'>mean</font>, "
    "<font color='#c0c0c0'>min</font>)");
  text.setFont(QFont("Arial", 10));
  plot_->setAxisTitle(QwtPlot::yRight, text);
}

void ErrorPlotWidget::setGenerationAxis() {
  QwtText text("Generation");
  text.setFont(QFont("Arial", 10));
  plot_->setAxisTitle(QwtPlot::xBottom, text);
}

void ErrorPlotWidget::setTimeAxis() {
  QwtText text("Time (seconds)");
  text.setFont(QFont("Arial", 10));
  plot_->setAxisTitle(QwtPlot::xBottom, text);
}

void ErrorPlotWidget::setLinScale() {
  plot_->setAxisScaleEngine(QwtPlot::xBottom, new QwtLinearScaleEngine());
  plot_->replot();
}

void ErrorPlotWidget::setLogScale() {
  plot_->setAxisScaleEngine(QwtPlot::xBottom, new QwtLogScaleEngine());
  plot_->replot();
}

void ErrorPlotWidget::clear() {
  plot_->detachItems();

  int n = valuesList.size();
  for (int i = 0; i < n; ++i)
    delete[] valuesList[i];

  valuesList.clear();

  plot_->replot();
}

// TODO: restrict to a single searchId
int ErrorPlotWidget::plotTimeBest(int searchId, DB *db) {
  clear();
  setSimpleErrorAxis();
  setSimpleCompAxis();
  setTimeAxis();
  int time = addPlotTime("", db);
  plot_->replot();
  return time;
}

int ErrorPlotWidget::plotTimeAll(int searchId, DB *db) {
  clear();
  setSimpleErrorAxis();
  setSimpleCompAxis();
  setTimeAxis();

  QSqlQuery *query = db->newTableQuery("Deme", "Search", searchId);

  int maxTime = 0;
  while (query->next()) {
    int demeId = query->value(0).toInt();
    int time = addPlotTime(QString("WHERE Generation.Deme = %1")
      .arg(demeId), db);
    if (time > maxTime)
      maxTime = time;
  }

  delete query;

  plot_->replot();
  return maxTime;
}

int ErrorPlotWidget::plotTimeDeme(int demeId, DB *db) {
  clear();
  setSimpleErrorAxis();
  setSimpleCompAxis();
  setTimeAxis();

  int time = addPlotTime(QString("WHERE Generation.Deme = %1").arg(demeId), db);

  plot_->replot();
  return time;
}

int ErrorPlotWidget::addPlotTime(const QString &cond, DB *db) {
  QSqlQuery *query = db->newQuery("SELECT MIN(Generation.Time), "
    "MAX(Generation.Time) FROM Generation " + cond);

  if (!query->next())
    return 0;

  int minTime = query->value(0).toInt();
  int maxTime = query->value(1).toInt();
  delete query;

  query = db->newQuery("SELECT Generation.Time, Generation.MinError, "
    "Generation.BestComp FROM Generation " + cond +
    " ORDER BY Generation.Time");

  double *errorXAxis = new double[maxTime - minTime + 1];
  valuesList.append(errorXAxis);
  double *compXAxis = new double[maxTime - minTime + 1];
  valuesList.append(compXAxis);
  double *errorYAxis = new double[maxTime - minTime + 1];
  valuesList.append(errorYAxis);
  double *compYAxis = new double[maxTime - minTime + 1];
  valuesList.append(compYAxis);

  double minError = 1e100;
  int minErrorTime = -1;
  int bestComp = 1e10;
  int bestCompTime = -1;
  int iError = -1;
  int iComp = -1;
  int time;
  while (query->next()) {
    time = query->value(0).toInt();
    double error = query->value(1).toDouble();
    int comp = query->value(2).toInt();

    // Accumulative min error and complexity at that error
    if ((error < minError && comp != bestComp) ||
      (error == minError && comp < bestComp)) {
      if (bestCompTime < time) { // New time point necessary
        iComp++;
        compXAxis[iComp] = time;
        bestCompTime = time;
      }

      compYAxis[iComp] = comp;
      bestComp = comp;
    }

    if (error < minError) {
      if (minErrorTime < time) { // New time point necessary
        iError++;
        errorXAxis[iError] = time;
        minErrorTime = time;
      }

      errorYAxis[iError] = error;
      minError = error;
    }
  }

  // Last points
  if (minErrorTime < time) {
    iError++;
    errorXAxis[iError] = time;
    errorYAxis[iError] = minError;
  }

  if (bestCompTime < time) {
    iComp++;
    compXAxis[iComp] = time;
    compYAxis[iComp] = bestComp;
  }

  delete query;

  QwtPlotCurve *curve;

  curve = new QwtPlotCurve();
  curve->setStyle(QwtPlotCurve::Steps);
  curve->setCurveAttribute(QwtPlotCurve::Inverted);
  curve->setPen(QPen(Qt::lightGray));
  //curve->setRenderHint(QwtPlotItem::RenderAntialiased);
  curve->setRawSamples(compXAxis, compYAxis, iComp + 1);
  curve->setYAxis(QwtPlot::yRight);
  curve->setZ(0);
  curve->attach(plot_);

  curve = new QwtPlotCurve();
  curve->setStyle(QwtPlotCurve::Steps);
  curve->setCurveAttribute(QwtPlotCurve::Inverted);
  curve->setPen(QPen(Qt::blue));
  //curve->setRenderHint(QwtPlotItem::RenderAntialiased);
  curve->setRawSamples(errorXAxis, errorYAxis, iError + 1);
  curve->setZ(1);
  curve->attach(plot_);

  return time;
}

// TODO: restrict to a single searchId
int ErrorPlotWidget::plotGeneBest(int searchId, DB *db) {
  clear();
  setFullErrorAxis();
  setFullCompAxis();
  setGenerationAxis();
  int maxGene = addPlotGene("", db);
  plot_->replot();
  return maxGene;
}

int ErrorPlotWidget::plotGeneAll(int searchId, DB *db) {
  clear();
  setFullErrorAxis();
  setFullCompAxis();
  setGenerationAxis();

  QSqlQuery *query = db->newTableQuery("Deme", "Search", searchId);

  int maxGene = 0;
  while (query->next()) {
    int demeId = query->value(0).toInt();
    int gene = addPlotGene(QString("WHERE Generation.Deme = %1")
      .arg(demeId), db);
    if (gene > maxGene)
      maxGene = gene;
  }

  delete query;

  plot_->replot();
  return maxGene;
}

int ErrorPlotWidget::plotGeneDeme(int demeId, DB *db) {
  clear();
  setFullErrorAxis();
  setFullCompAxis();
  setGenerationAxis();
  int maxGene = addPlotGene(QString("WHERE Generation.Deme = %1")
    .arg(demeId), db);

  plot_->replot();
  return maxGene;
}

int ErrorPlotWidget::addPlotGene(const QString &cond, DB *db) {
  QSqlQuery *query;

  if (cond.isEmpty())
    query = db->newQuery("SELECT 1 + MAX(Generation.Ind) FROM Generation");
  else
    query = db->newQuery("SELECT COUNT(Generation.Ind) FROM Generation " + cond);

  if (!query->next())
    return 0;

  int nGene = query->value(0).toInt();
  delete query;

  if (cond.isEmpty()) // Best deme
    query = db->newQuery("SELECT Generation.MinError, "
      "AVG(Generation.MeanError), MAX(Generation.MaxError), "
      "MIN(Generation.MinComp), AVG(Generation.MeanComp), "
      "MAX(Generation.MaxComp), MIN(Generation.BestComp) "
      "FROM Generation INNER JOIN "
      "(SELECT Generation.Ind, MIN(Generation.MinError) AS MinError "
      "FROM Generation GROUP BY Generation.Ind) AS sub_query "
      "ON Generation.Ind = sub_query.Ind AND "
      "Generation.MinError = sub_query.MinError GROUP BY Generation.Ind "
      "ORDER BY Generation.Ind");
  else
    query = db->newQuery("SELECT Generation.MinError, "
      "Generation.MeanError, Generation.MaxError, Generation.MinComp, "
      "Generation.MeanComp, Generation.MaxComp, Generation.BestComp "
      "FROM Generation " + cond);

  double *xAxis = new double[nGene];
  valuesList.append(xAxis);

  for (int i = 0; i < nGene; ++i)
    xAxis[i] = i + 1;


  double *minError = new double[nGene];
  valuesList.append(minError);
  double *meanError = new double[nGene];
  valuesList.append(meanError);
  double *maxError = new double[nGene];
  valuesList.append(maxError);
  double *minComp = new double[nGene];
  valuesList.append(minComp);
  double *meanComp = new double[nGene];
  valuesList.append(meanComp);
  double *maxComp = new double[nGene];
  valuesList.append(maxComp);
  double *bestComp = new double[nGene];
  valuesList.append(bestComp);

  int i = 0;
  while (query->next()) {
    minError[i] = query->value(0).toDouble();
    meanError[i] = query->value(1).toDouble();
    maxError[i] = query->value(2).toDouble();
    minComp[i] = query->value(3).toDouble();
    meanComp[i] = query->value(4).toDouble();
    maxComp[i] = query->value(5).toDouble();
    bestComp[i] = query->value(6).toDouble();
    ++i;
  }

  delete query;

  QwtPlotCurve *curve;

  curve = new QwtPlotCurve();
  curve->setStyle(QwtPlotCurve::Steps);
  curve->setCurveAttribute(QwtPlotCurve::Inverted);
  curve->setPen(QPen(Qt::lightGray));
  //curve->setRenderHint(QwtPlotItem::RenderAntialiased);
  curve->setRawSamples(xAxis, minComp, nGene);
  curve->setYAxis(QwtPlot::yRight);
  curve->setZ(0);
  curve->attach(plot_);

  curve = new QwtPlotCurve();
  curve->setStyle(QwtPlotCurve::Steps);
  curve->setCurveAttribute(QwtPlotCurve::Inverted);
  curve->setPen(QPen(Qt::lightGray));
  //curve->setRenderHint(QwtPlotItem::RenderAntialiased);
  curve->setRawSamples(xAxis, meanComp, nGene);
  curve->setYAxis(QwtPlot::yRight);
  curve->setZ(0);
  curve->attach(plot_);

  curve = new QwtPlotCurve();
  curve->setStyle(QwtPlotCurve::Steps);
  curve->setCurveAttribute(QwtPlotCurve::Inverted);
  curve->setPen(QPen(Qt::lightGray));
  //curve->setRenderHint(QwtPlotItem::RenderAntialiased);
  curve->setRawSamples(xAxis, maxComp, nGene);
  curve->setYAxis(QwtPlot::yRight);
  curve->setZ(0);
  curve->attach(plot_);

  curve = new QwtPlotCurve();
  curve->setStyle(QwtPlotCurve::Steps);
  curve->setCurveAttribute(QwtPlotCurve::Inverted);
  curve->setPen(QPen(Qt::black));
  //curve->setRenderHint(QwtPlotItem::RenderAntialiased);
  curve->setRawSamples(xAxis, bestComp, nGene);
  curve->setYAxis(QwtPlot::yRight);
  curve->setZ(1);
  curve->attach(plot_);

  curve = new QwtPlotCurve();
  curve->setStyle(QwtPlotCurve::Steps);
  curve->setCurveAttribute(QwtPlotCurve::Inverted);
  curve->setPen(QPen(Qt::red));
  //curve->setRenderHint(QwtPlotItem::RenderAntialiased);
  curve->setRawSamples(xAxis, maxError, nGene);
  curve->setZ(2);
  curve->attach(plot_);

  curve = new QwtPlotCurve();
  curve->setStyle(QwtPlotCurve::Steps);
  curve->setCurveAttribute(QwtPlotCurve::Inverted);
  curve->setPen(QPen(Qt::green));
  //curve->setRenderHint(QwtPlotItem::RenderAntialiased);
  curve->setRawSamples(xAxis, meanError, nGene);
  curve->setZ(3);
  curve->attach(plot_);

  curve = new QwtPlotCurve();
  curve->setStyle(QwtPlotCurve::Steps);
  curve->setCurveAttribute(QwtPlotCurve::Inverted);
  curve->setPen(QPen(Qt::blue));
  //curve->setRenderHint(QwtPlotItem::RenderAntialiased);
  curve->setRawSamples(xAxis, minError, nGene);
  curve->setZ(4);
  curve->attach(plot_);


  return nGene;
}

void ErrorPlotWidget::PlotMagnifierBottomFixed::rescale(double factor) {
  factor = qAbs(factor);
  if (factor == 1.0 || factor == 0.0)
    return;

  bool doReplot = false;
  QwtPlot* plt = plot();

  const bool autoReplot = plt->autoReplot();
  plt->setAutoReplot(false);

  for (int axisId = 0; axisId < QwtPlot::axisCnt; axisId++) {
    const QwtScaleDiv &scaleDiv = plt->axisScaleDiv(axisId);
    if (isAxisEnabled(axisId) && !scaleDiv.isEmpty()) {
      double lb = scaleDiv.lowerBound();
      double ub = scaleDiv.upperBound();
      plt->setAxisScale(axisId, lb, lb + (ub - lb) * factor);
      doReplot = true;
    }
  }

  plt->setAutoReplot(autoReplot);

  if (doReplot)
    plt->replot();
}

QwtText ErrorPlotWidget::Power10ScaleDraw::label(double v) const {
  return QString("10<sup>%1</sup>").arg(log10(v));
}

// private slot
void ErrorPlotWidget::saveImage() {
  QSettings settings;
  settings.beginGroup("ErrorPlotWidget");
  QString lastDir = settings.value("lastDir").toString();

  QString fileName = QFileDialog::getSaveFileName(this, "Save image",
    lastDir, "Images (*.svg *.pdf *.ps *.png *.jpg *.bmp)");

  if (!fileName.isEmpty()) {
    QString newDir = QFileInfo(fileName).absolutePath();
    if (newDir != lastDir)
      settings.setValue("lastDir", newDir);

    QwtPlotRenderer renderer;
    //plot_->replot();

    QSizeF renderSize = size();
    renderSize.scale(100, 100, Qt::KeepAspectRatio);
    renderer.renderDocument(plot_, fileName, renderSize, 100);
  }
}

}