// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#pragma once

#include "Private/zoomgraphicsview.h"

#include <QGraphicsSvgItem>
#include <QProcess>
#include <QTextStream>

namespace LoboLab {

class Model;
class ModelComp;

class ModelGraphView : public ZoomGraphicsView {
  Q_OBJECT

  public:
    ModelGraphView(Model *m, int nMorphogens, bool hideNotUsed = false, QWidget * parent = NULL);
    virtual ~ModelGraphView();

    public slots:
    void updateModel(bool hideNotUsed);

  signals:
    //void modified();

  private:
    void createGraph(bool hideNotUsed);
    void writeModel2Dot(QTextStream &stream, bool hideNotUsed);
    void writeProd(QTextStream &stream, int label, int fillColorInd,
      const QString name) const;
    void writeProdAttribs(QTextStream &stream, const QStringList &props,
      const QStringList &style) const;

    Model *model_;
    QGraphicsSvgItem *svgItem_;
    QProcess *process_;
    QStringList args_;

    int nMorphogens_;
};

} // namespace LoboLab
