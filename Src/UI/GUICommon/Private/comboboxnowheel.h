// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#pragma once

#include <QComboBox>

namespace LoboLab {

class ComboBoxNoWheel : public QComboBox {
 public:
  ComboBoxNoWheel(QWidget * parent = 0);
  virtual ~ComboBoxNoWheel();

  virtual bool event(QEvent *event);
};

} // namespace LoboLab
