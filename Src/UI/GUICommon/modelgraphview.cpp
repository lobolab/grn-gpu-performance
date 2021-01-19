// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#include "modelgraphview.h"

#include "Private/colorsconfig.h"
#include "Model/model.h"
#include "Model/modelProd.h"
#include "Model/modelLink.h"

#include <QFile>
#include <QFileInfo>
#include <QStringList>

namespace LoboLab {

ModelGraphView::ModelGraphView(Model *m, int nMorphogens, bool hideNotUsed, QWidget * parent)
  : ZoomGraphicsView(parent, true, false),
  model_(m),
  nMorphogens_(nMorphogens) {
  process_ = new QProcess(this);

  args_.append("-Tsvg");
  args_.append("graph.gv");
  args_.append("-ograph.svg");

  createGraph(hideNotUsed);
}

ModelGraphView::~ModelGraphView() {
  delete process_;
}

void ModelGraphView::createGraph(bool hideNotUsed) {
  QFile gvFile("graph.gv");

  gvFile.open(QIODevice::Truncate | QIODevice::WriteOnly | QIODevice::Text);
  Q_ASSERT(gvFile.isOpen());

  QTextStream gvStream(&gvFile);
  writeModel2Dot(gvStream, hideNotUsed);
  gvStream.flush();
  gvFile.close();

  process_->start("dot", args_);
  process_->waitForFinished();

  QFile svgFile("graph.svg");
  Q_ASSERT(svgFile.exists());

  svgItem_ = new QGraphicsSvgItem("graph.svg");

  scene_->addItem(svgItem_);

  bool ok = gvFile.remove();
  Q_ASSERT(ok);

  ok = QFile::remove("graph.svg");
  Q_ASSERT(ok);
}

void ModelGraphView::writeModel2Dot(QTextStream &stream, bool hideNotUsed) {
  stream << "digraph G {bgcolor=transparent;";
  stream << "node [fontname=Arial];";

  QList<int> labels = model_->calcProductLabels().toList();
  qSort(labels);
  // Part of GUI. did not implement to take Number of Morphogens from the database.
  QSet<int> usedLabels = model_->calcProductLabelsInUse(nMorphogens_);
  int nProducts = labels.size();

  QStringList names;
  if (nMorphogens_ >= 1)
    names.append("Red");
  if (nMorphogens_ >= 2)
    names.append("Green");
  if (nMorphogens_ >= 3)
    names.append("\" Blue  \"");
  
  int c = 97;
  for (int i = names.size(); i < nProducts; ++i) {
    if (usedLabels.contains(labels.at(i))) {
      names.append(QString("%1"/*"<<i>%1</i>>"*/).arg(QChar(c)));
      ++c;
    }
    else {
      names.append(QString("%1").arg(labels.at(i)));
    }
  }

  int colorInd = 0;
  for (int i = 0; i < nProducts; ++i) {
    if (usedLabels.contains(labels.at(i)))
      writeProd(stream, labels.at(i), colorInd++, names[i]);
    else if (!hideNotUsed)
      writeProd(stream, labels.at(i), -1, names[i]);
  }

  // Links
  int n = labels.size();
  for (int i = 0; i < n; ++i) {
    int regulatedLabel = labels.at(i);
    QList<ModelLink*> links = model_->linksToLabel(regulatedLabel);

    int nLinks = links.size();
    for (int j = 0; j < nLinks; ++j) {
      ModelLink *link = links.at(j);

      bool linkIsUsed = usedLabels.contains(link->regulatorProdLabel()) &&
        usedLabels.contains(link->regulatedProdLabel());

      if (linkIsUsed || !hideNotUsed) {
        stream << link->regulatorProdLabel() << " -> " << regulatedLabel;

        QString color;
        QString arrowhead;

        if (link->hillCoef() < 0) {
          color = "ff0000";
          arrowhead = "tee";
        }
        else {
          color = "0000ff";
          arrowhead = "normal";
        }

        if (!linkIsUsed)
          color += "80"; // alpha channel

        stream << " [color=\"#" << color << "\",arrowhead=" << arrowhead;

        if (link->isAndReg())
          stream << ",style=bold];";
        else
          stream << ",style=\"bold,dashed\"];";
      }
    }
  }

  stream << "}";
}

void ModelGraphView::writeProd(QTextStream &stream, int label, int fillColorInd,
  const QString name) const {
  stream << label;

  QStringList props;
  QStringList style;

  props += QString("label=%1").arg(name);

  if (label < nMorphogens_) {
    props += "shape=diamond";
    props += "margin=\"0.05,0.076\"";
  }

  if (fillColorInd > -1) {
    style += "filled";
    // With alpha channel
    QRgb color = 0x80 + 0x100 * (
      ColorsConfig::colorsCells[fillColorInd % ColorsConfig::nColors]);
    // Whithout alpha channel
    //QRgb color = 0xFF + 0x100 * (
    //  ColorsConfig::colorsCells[fillColorInd % ColorsConfig::nColors]);
    QString colStr = QString("\"#%1\"").arg(color, 8, 16, QChar('0')).toUpper();
    props += "fillcolor=" + colStr;

    if (model_->prodWithLabel(label)->dif() > 0)
      //props += "color=" + colStr;
      style += "dashed";
    //else
    //props += "color=black";
  }
  else {
    props += "fontcolor=gray";
    props += "color=gray";

    if (model_->prodWithLabel(label)->dif() > 0)
      style += "dashed";
  }

  writeProdAttribs(stream, props, style);

  stream << ';';
}

void ModelGraphView::writeProdAttribs(QTextStream &stream,
  const QStringList &props, const QStringList &style) const {
  if (!props.isEmpty() || !style.isEmpty()) {
    stream << " [";

    if (!props.isEmpty()) {
      stream << props.first();
      for (int i = 1; i < props.size(); ++i)
        stream << ',' << props.at(i);
    }

    if (!style.isEmpty()) {
      if (!props.isEmpty())
        stream << ',';

      stream << "style=\"";

      stream << style.first();
      for (int i = 1; i < style.size(); ++i)
        stream << ',' << style.at(i);

      stream << '\"';
    }

    stream << ']';
  }
}

void ModelGraphView::updateModel(bool hideNotUsed) {
  delete svgItem_;
  createGraph(hideNotUsed);
  fitView();
}


}