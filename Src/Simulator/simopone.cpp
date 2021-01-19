// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#include "simopone.h"

namespace LoboLab {

SimOpOne::SimOpOne(double *t)
  : to_(t) {
}

SimOpOne::~SimOpOne() {
}

void SimOpOne::compute() const {
  *to_ = 1;
}

}