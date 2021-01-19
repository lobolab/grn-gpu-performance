// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#pragma once

namespace LoboLab {
class CSimOp {
  public: 
	  typedef enum { OpAnd, OpOr, OpHillAct, OpHillRep, OpZero } FuncType;
	  typedef union { int index;  double *pointer; } IndexOrPointerToData;

		__host__ __device__ CSimOp();
		__host__ __device__ CSimOp(FuncType funcType, int to, int from = 0, double dc = 0.0, double hc = 0.0);
		__device__ void compute();

		template <class T> __device__ static inline const T& min(const T& a, const T& b);
		__device__ static void and_func(CSimOp *op);

		template <class T> __device__ static inline const T& max(const T& a, const T& b);
		__device__ static void or_func(CSimOp *op);

		__device__ static void zero_func(CSimOp *op);
		__device__ static void hillAct_func(CSimOp *op);
		__device__ static void hillRep_func(CSimOp *op);
		__device__ void linkFuncPointer(double *ratios, double *oldConc, double *tempHill, int nThreads);
     
  private:

    //See https://devtalk.nvidia.com/default/topic/821546/math-pow-double-2-0-has-extremely-low-branches-efficiency-but-math-pow-double-2-is-fast-/
    __host__ __device__ inline static double powTrick(double x, double y) {
      return exp(y*log(x));
    }

    void(*op_func_) (CSimOp*);
		
	  FuncType funcType_;

	  IndexOrPointerToData from_;
	  IndexOrPointerToData to_;

	  double hillCoef_;
	  double disConstN_;
};

}