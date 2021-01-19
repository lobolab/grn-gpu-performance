// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#include "simophillrep.h"
#include "Common/mathalgo.h"

namespace LoboLab {

SimOpHillRep::SimOpHillRep(const double *f, double *t, double dc, double hc)
  : from_(f), to_(t), disConst_(dc), hillCoef_(hc),
    disConstN_(pow(disConst_, hillCoef_)) {
}

SimOpHillRep::~SimOpHillRep() {
}

void SimOpHillRep::compute() const {
  double rcn = pow(*from_, hillCoef_);
  *to_ = disConstN_ / (rcn + disConstN_);
}

}