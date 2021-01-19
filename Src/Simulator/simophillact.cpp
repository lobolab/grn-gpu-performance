// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#include "simophillact.h"
#include "Common/mathalgo.h"

namespace LoboLab {

SimOpHillAct::SimOpHillAct(const double *f, double *t, double dc, double hc)
  : from_(f), to_(t), disConst_(dc), hillCoef_(hc),
    disConstN_(pow(disConst_, hillCoef_)) {
}

SimOpHillAct::~SimOpHillAct() {
}

void SimOpHillAct::compute() const {
  double rcn = pow(*from_, hillCoef_);
  *to_ = rcn / (rcn + disConstN_);
}

}