// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#pragma once

#include "simop.h"

namespace LoboLab {

class SimOpOne : public SimOp {
 public:
  explicit SimOpOne(double *to);
  ~SimOpOne();

  void compute() const;

 private:
  double *to_;
};

} // namespace LoboLab
