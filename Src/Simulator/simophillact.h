// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#pragma once

#include "simop.h"

namespace LoboLab {

class SimOpHillAct : public SimOp {
 public:
  SimOpHillAct(const double *from, double *to, double disConst, double hillCoef);
  ~SimOpHillAct();

  void compute() const;

 private:
  const double *from_;
  double *to_;
  const double disConst_;
  const double hillCoef_;
  const double disConstN_;
};

} // namespace LoboLab
