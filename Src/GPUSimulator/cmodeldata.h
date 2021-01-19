// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#pragma once

namespace LoboLab {

class CSimOp;
struct CSimState;

struct CModelData {

  int nProducts;

  CSimOp* ops;
  int nOps;
  double *limits;
  double *degradations;
  int *difProdInd;
  double *difConsts;
  int nDif;

  double maxError;
};

}