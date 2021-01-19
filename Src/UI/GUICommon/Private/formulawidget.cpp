// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#include "formulawidget.h"

#include <QPaintEvent>
#include <QDebug>

#include <qwt_mml_document.h>

namespace LoboLab {

FormulaWidget::FormulaWidget(QWidget *parent)
  : QLabel(parent),
    backgroundColor_(Qt::white),
    foregroundColor_(Qt::black),
    fontSize_(11),
    rotation_(0) {
  setStyleSheet("background-color:white;");
  setAlignment(Qt::AlignLeft | Qt::AlignTop);
  setMargin(5);
}

FormulaWidget::~FormulaWidget() {
}

void FormulaWidget::setMathMLStr(const QString &formulaStr) {
  formulaStr_ = formulaStr;
  paintFormula();
  update();
}

void FormulaWidget::setFontSize(double fontSize) {
  fontSize_ = fontSize;
  update();
}

void FormulaWidget::setRotation(double rotation) {
  rotation_ = rotation;
  update();
}

void FormulaWidget::setBackgroundColor(const QColor &color) {
  backgroundColor_ = color;
  update();
}

void FormulaWidget::setForegroundColor(const QColor &color) {
  foregroundColor_ = color;
  update();
}

void FormulaWidget::paintFormula() {  
  QwtMathMLDocument doc;
  doc.setContent(formulaStr_);

  doc.setBaseFontPointSize(fontSize_);
  
  QSize docSize = doc.size();

  
  QPixmap pixmap(docSize);
  pixmap.fill(Qt::white);
  
  QPainter painter(&pixmap);
  
  doc.paint(&painter, QPoint(0, 0));

  painter.end();

  setPixmap(pixmap);
  setMinimumSize(100,100);
}

}
