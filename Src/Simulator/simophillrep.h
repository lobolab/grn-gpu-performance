// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#pragma once

#include "simop.h"

namespace LoboLab {

class SimOpHillRep : public SimOp {
 public:
  SimOpHillRep(const double *from, double *to, double disConst, double hillCoef);
  ~SimOpHillRep();

  void compute() const;

 private:
  const double *from_;
  double *to_;
  double disConst_;
  double hillCoef_;
  double disConstN_;
};

} // namespace LoboLab
