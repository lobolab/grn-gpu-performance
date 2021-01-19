// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#include "doublespinboxnowheel.h"

#include <QEvent>

namespace LoboLab {

DoubleSpinBoxNoWheel::DoubleSpinBoxNoWheel(QWidget * parent)
  : QDoubleSpinBox(parent) {
}


DoubleSpinBoxNoWheel::~DoubleSpinBoxNoWheel() {
}

bool DoubleSpinBoxNoWheel::event(QEvent * event) {
  if (event->type() == QEvent::Wheel)
    return false;
  else
    return QDoubleSpinBox::event(event);
}

}