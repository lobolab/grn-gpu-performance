// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#pragma once

#include "simop.h"

namespace LoboLab {

class SimOpOr : public SimOp {
 public:
  SimOpOr(const double *from1, double *from2AndTo);
  SimOpOr(const double *from1, const double *from2, double *to);
  ~SimOpOr();

  void compute() const;

 private:
  const double *from1_;
  const double *from2_;
  double *to_;
};

} // namespace LoboLab
