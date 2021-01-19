// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#pragma once

namespace LoboLab {

struct CSimState;

struct CExperimentsData {
    
  int nExperiments;
  int nMaxSteps;
  double minConc;
  double minConcChange;
  double dt;
  double distErrorThreshold;
  int nInputMorphogens;
  int nTargetMorphogens;

  CSimState ** inSimStates;
  CSimState ** targetSimStates;
};

}