// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#include "comboboxnowheel.h"

#include <QEvent>

namespace LoboLab {

ComboBoxNoWheel::ComboBoxNoWheel(QWidget * parent)
  : QComboBox(parent) {
}


ComboBoxNoWheel::~ComboBoxNoWheel() {
}

bool ComboBoxNoWheel::event(QEvent * event) {
  if (event->type() == QEvent::Wheel)
    return false;
  else
    return QComboBox::event(event);
}

}