// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#include "simopand.h"
#include "Common/mathalgo.h"

namespace LoboLab {

SimOpAnd::SimOpAnd(const double *f1, double *f2AndT)
  : from1_(f1), from2_(f2AndT), to_(f2AndT) {
}

SimOpAnd::SimOpAnd(const double *f1, const double *f2, double *t)
  : from1_(f1), from2_(f2), to_(t) {
}

SimOpAnd::~SimOpAnd() {
}

void SimOpAnd::compute() const {
  *to_ = MathAlgo::min(*from1_, *from2_);
}

}