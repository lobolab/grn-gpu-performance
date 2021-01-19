// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#include "simopcopy.h"

namespace LoboLab {

SimOpCopy::SimOpCopy(const double *f, double *t)
  : from_(f), to_(t) {
}

SimOpCopy::~SimOpCopy() {
}

void SimOpCopy::compute() const {
  *to_ = *from_;
}

}