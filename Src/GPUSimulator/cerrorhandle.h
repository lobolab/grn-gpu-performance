// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#pragma once

#include "Common/log.h"

#include <stdio.h>
#include <stdlib.h>

#include <QString>

#include <cuda_runtime.h>
#include <device_launch_parameters.h>

namespace LoboLab {

#ifdef QT_DEBUG
#define cErrorHandle(gpu_function_return) {cudaError_t gpu_function_return_code = (gpu_function_return); Q_ASSERT_X(gpu_function_return_code == cudaSuccess, ("code error " + QString::number(gpu_function_return_code)).toLatin1().data(), cudaGetErrorString(gpu_function_return_code));}
#else
#define cErrorHandle(gpu_function_return) {cudaError_t gpu_function_return_code = (gpu_function_return); if(gpu_function_return_code != cudaSuccess) {Log::write() << "GPU ERROR: " << "Code: " << QString::number(gpu_function_return_code) << " Msg: " << cudaGetErrorString(gpu_function_return_code) << endl; abort();}}
#endif

#ifdef QT_DEBUG
#define GPU_ASSERT(condition) if (!(condition)) { printf("ERROR: GPU_ASSERT failed: tid %d: %s, %d\n", threadIdx.x, __FILE__, __LINE__); asm("trap;"); }
#else
#define GPU_ASSERT( x )
#endif
}