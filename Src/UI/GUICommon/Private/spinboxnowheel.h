// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#pragma once

#include <QSpinBox>

namespace LoboLab {

class SpinBoxNoWheel : public QSpinBox {
 public:
  SpinBoxNoWheel(QWidget * parent = 0);
  virtual ~SpinBoxNoWheel();

  virtual bool event(QEvent *event);
};

} // namespace LoboLab
