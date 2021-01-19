// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#pragma once

#include <cuda_runtime.h>

namespace LoboLab {

  struct CSimState {

	  __host__ void initialize(int w, int h, int nProducts);
	  __host__ void freeMatrix();
	  __device__ double getProduct(int row, int column, int productIndx);
	  __device__ void setProduct(int row, int column, int productIndx, double productVal);
	  __host__ void resize(int newNumProducts);
	  __host__ CSimState* copyToDevice(cudaStream_t stream);

	  __device__ __host__ inline int nProducts() { return nProducts_; }
	  __device__ __host__ inline int width() { return width_; }
	  __device__ __host__ inline int height(){ return height_; }
	  
	  int nProducts_;
	  int width_;
	  int height_;
	  double *matrix_;
  };
}

    