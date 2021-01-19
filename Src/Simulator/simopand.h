// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#pragma once

#include "simop.h"

namespace LoboLab {

class SimOpAnd : public SimOp {
 public:
  SimOpAnd(const double *from1, double *from2AndTo);
  SimOpAnd(const double *from1, const double *from2, double *to);
  ~SimOpAnd();

  void compute() const;

 private:
  const double *from1_;
  const double *from2_;
  double *to_;
};

} // namespace LoboLab
