// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#pragma once

namespace LoboLab {

struct CSimTempData {
  CSimState cSimState1;
  CSimState cSimState2;

  double *ratios;
  double *oldConcs;
};

}