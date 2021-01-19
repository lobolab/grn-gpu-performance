// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#pragma once

#include "simop.h"

namespace LoboLab {

class SimOpCopy : public SimOp {
 public:
  SimOpCopy(const double *from, double *to);
  ~SimOpCopy();

  void compute() const;

 private:
  const double *from_;
  double *to_;
};

} // namespace LoboLab
