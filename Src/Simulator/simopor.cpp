// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#include "simopor.h"
#include "Common/mathalgo.h"

namespace LoboLab {

SimOpOr::SimOpOr(const double *f1, double *f2AndT)
  : from1_(f1), from2_(f2AndT), to_(f2AndT) {
}

SimOpOr::SimOpOr(const double *f1, const double *f2, double *t)
  : from1_(f1), from2_(f2), to_(t) {
}

SimOpOr::~SimOpOr() {
}

void SimOpOr::compute() const {
  *to_ = MathAlgo::max(*from1_, *from2_);
}

}