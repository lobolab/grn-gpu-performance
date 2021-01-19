// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#include "csimstate.cuh"
#include "cerrorhandle.h"

namespace LoboLab {

__host__ void CSimState::initialize(int width, int height, int nProducts) {
  width_ = width;
  height_ = height;
  nProducts_ = nProducts;
		
  int totalSize = width_ * height_ * nProducts_;

  if (nProducts_ > 0) {
    cErrorHandle(cudaMalloc(&matrix_, sizeof(double)*(totalSize)));
  }
  else {
    matrix_ = NULL;
  }
}

__host__ void CSimState::freeMatrix() {
	if (matrix_)
    cErrorHandle(cudaFree(matrix_));

  matrix_ = NULL;
}

__device__ double CSimState::getProduct(int row, int column, int productIndx) {
	int index = productIndx * width_ * height_ + row * width_ + column;
  GPU_ASSERT(index >= 0 && index < (width_ * height_ * nProducts_));
  return matrix_[index];
}

__device__ void CSimState::setProduct(int row, int column, int productIndx, double productVal) {
	int index = productIndx * width_ * height_ + row * width_ + column;
  GPU_ASSERT(index >= 0 && index < (width_ * height_ * nProducts_));
	matrix_[index] = productVal;
}

__host__ void CSimState::resize(int newNProducts) {
	if (newNProducts != nProducts_) {
    cErrorHandle(cudaFree(matrix_));
    nProducts_ = newNProducts;
		int totalSize_ = width_ * height_ * nProducts_;
    cErrorHandle(cudaMalloc(&matrix_, sizeof(double)*(totalSize_)));
	}
}
}