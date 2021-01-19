// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#pragma once

#include "simop.h"

namespace LoboLab {

class SimOpZero : public SimOp {
 public:
  explicit SimOpZero(double *to);
  ~SimOpZero();

  void compute() const;

 private:
  double *to_;
};

} // namespace LoboLab
