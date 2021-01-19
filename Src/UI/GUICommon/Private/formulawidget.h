// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.
// Inspired by https://github.com/uwerat/qwt-mml-dev

#pragma once

#include <QLabel>
#include <QPainter>

namespace LoboLab {

class FormulaWidget: public QLabel {
  Q_OBJECT

  public:
    FormulaWidget(QWidget *parent = NULL);
    virtual ~FormulaWidget();

    inline const QString &mathMLStr() const { return formulaStr_; };

    void setMathMLStr(const QString & formulaStr);
    void setBackgroundColor(const QColor &color);
    void setForegroundColor(const QColor &color);
    void setFontSize(double fontSize);
    void setRotation(double rotation);
    
  private:
    void paintFormula();

  private:
    QString formulaStr_;
    QColor backgroundColor_;
    QColor foregroundColor_;
    double fontSize_;
    double rotation_;
};

}
