// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#include "csimop.cuh"

#include <stdio.h>

//#include "device_launch_parameters.h"

namespace LoboLab {

  __host__ __device__ CSimOp::CSimOp() {}

  __host__ __device__ CSimOp::CSimOp(FuncType funcType, int to, int from, double dc, double hc)
    : op_func_(NULL), funcType_(funcType), from_{from}, to_{to}, hillCoef_(hc), disConstN_(CSimOp::powTrick(dc, hc)) {}

  __device__ void CSimOp::compute() {
      op_func_(this);
  }

  template <class T>
  __device__ inline const T& CSimOp::min(const T& a, const T& b) {
    return (!(b < a)) ? a : b;
  }

  __device__ void CSimOp::and_func(CSimOp *op) {
    op->to_.pointer[threadIdx.x] = min(op->from_.pointer[threadIdx.x], op->to_.pointer[threadIdx.x]);
  }

  template <class T>
  __device__ inline const T& CSimOp::max(const T& a, const T& b) {
    return (a < b) ? b : a;
  }

  __device__ void CSimOp::or_func(CSimOp *op) {
    op->to_.pointer[threadIdx.x] = max(op->from_.pointer[threadIdx.x], op->to_.pointer[threadIdx.x]);
  }

  __device__ void CSimOp::zero_func(CSimOp *op) {
    op->to_.pointer[threadIdx.x] = 0.0;
  }
    
  __device__ void CSimOp::hillAct_func(CSimOp *op) {
    double rcn = CSimOp::powTrick(op->from_.pointer[threadIdx.x], op->hillCoef_);
    op->to_.pointer[threadIdx.x] = rcn / (rcn + op->disConstN_);
  }

  __device__ void CSimOp::hillRep_func(CSimOp *op) {
    double rcn = CSimOp::powTrick(op->from_.pointer[threadIdx.x], op->hillCoef_);
    op->to_.pointer[threadIdx.x] = op->disConstN_ / (rcn + op->disConstN_);
  }

  __device__ void CSimOp::linkFuncPointer(double *ratios, double *oldConc, double *tempHill, int nThreads) {

    switch (funcType_) {
    case OpAnd:
      op_func_ = and_func;
      from_.pointer = tempHill;
      to_.pointer = &ratios[to_.index*nThreads];
      break;
    case OpOr:
      op_func_ = or_func;
      from_.pointer = tempHill;
      to_.pointer = &ratios[to_.index*nThreads];
      break;
    case OpHillAct:
      op_func_ = hillAct_func;
      from_.pointer = &oldConc[from_.index*nThreads];
      if ((to_.index) < 0)
        to_.pointer = tempHill;
      else
        to_.pointer = &ratios[to_.index*nThreads];
      break;
    case OpHillRep:
      op_func_ = hillRep_func;
      from_.pointer = &oldConc[from_.index*nThreads];
      if ((to_.index) < 0)
        to_.pointer = tempHill;
      else
        to_.pointer = &ratios[to_.index*nThreads];
      break;
    case OpZero:
      op_func_ = zero_func;
      to_.pointer = &ratios[to_.index*nThreads];
      break;
    }
  }
}
