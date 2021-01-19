// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#pragma once

#include <QDoubleSpinBox>

namespace LoboLab {

class DoubleSpinBoxNoWheel : public QDoubleSpinBox {
 public:
  DoubleSpinBoxNoWheel(QWidget * parent = 0);
  virtual ~DoubleSpinBoxNoWheel();

  virtual bool event(QEvent *event);
};

} // namespace LoboLab
