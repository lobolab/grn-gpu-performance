// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#include "spinboxnowheel.h"

#include <QEvent>

namespace LoboLab {

SpinBoxNoWheel::SpinBoxNoWheel(QWidget * parent)
  : QSpinBox(parent) {
}

SpinBoxNoWheel::~SpinBoxNoWheel() {
}

bool SpinBoxNoWheel::event(QEvent * event) {
  if (event->type() == QEvent::Wheel)
    return false;
  else
    return QSpinBox::event(event);
}

}