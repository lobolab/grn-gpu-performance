// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#include "csimstate.cuh"
#include "modelsimulatordevice.cuh"
#include "cerrorhandle.h"
#include "cexperimentsdata.h"
#include "cmodeldata.h"
#include "csimtempdata.h"
#include "csimop.cuh"

#include "Common/log.h"

#include <QElapsedTimer>

#include <device_launch_parameters.h>
#include <device_functions.h>

namespace LoboLab {

// Storage in shared memory
struct CSharedData {
  double *tempHill;
  double *reductionArray;

  __device__ void initialize() {
    extern __shared__ double shared_mem[];
    tempHill = shared_mem; // nThreads 
    reductionArray = &shared_mem[blockDim.x]; // nThreads
  }
};


// Storage in local thread memory
struct CThreadData {
  int startIndex;
  int endIndex;
  int dimX;
  int dimY;

  __device__ void initialize(CExperimentsData *cExpData) {
    dimX = cExpData->inSimStates[0]->width();
    dimY = cExpData->inSimStates[0]->height();
    int nCellsPerThread = dimX * dimY / blockDim.x;
    int col = threadIdx.x;
    startIndex = col * nCellsPerThread;
    endIndex = startIndex + nCellsPerThread;
  }

};


__device__ void computeOps(CSimOp* ops, int nOps) {
  for (int i = 0; i < nOps; ++i)
    ops[i].compute();   
}


__device__ void computeRatios(CModelData *cModelData, CThreadData *cThreadData, 
    CSimTempData *cSimTempData, int row, int col, CSimState *simState) {
  for (int i = 0; i < cModelData->nProducts; ++i)
    cSimTempData->ratios[blockDim.x*i + threadIdx.x] = 
      cModelData->limits[i] * cSimTempData->ratios[blockDim.x*i + threadIdx.x] 
      - cModelData->degradations[i] * cSimTempData->oldConcs[blockDim.x*i + threadIdx.x];
    
  for (int d = 0; d < cModelData->nDif; ++d){
    double diffusion = 0.0;
    double total = 0.0;
    int k = cModelData->difProdInd[d];

    if (row > 0){
      diffusion += simState->getProduct(row - 1, col,k);
      ++total;
    }
    if (row < cThreadData->dimX - 1){
      diffusion += simState->getProduct(row + 1, col, k);
      ++total;
    }
    if (col > 0){
      diffusion += simState->getProduct(row, col - 1, k);
      ++total;
    }
    if (col < cThreadData->dimY - 1){
      diffusion += simState->getProduct(row, col + 1, k);
      ++total;
    }
      
    if (diffusion > 0 || (total > 0 && cSimTempData->oldConcs[blockDim.x*k + threadIdx.x] > 0))
      cSimTempData->ratios[blockDim.x*k + threadIdx.x] += 
        cModelData->difConsts[d] * 
        (diffusion - total * 
          cSimTempData->oldConcs[blockDim.x*k + threadIdx.x]);
  }
}


__device__ double computeError(CExperimentsData *cExpData, CSharedData *cSharedData, 
    CThreadData *cThreadData, CSimState* simState, CSimState* targetSimState) {
  double dist = 0.0;
  for (int i = cThreadData->startIndex; i < cThreadData->endIndex; ++i) {
    int row = i / cThreadData->dimX;
    int col = i % cThreadData->dimX;
    for (int k = 0; k < cExpData->nTargetMorphogens; ++k){
      double val1 = simState->getProduct(row, col, k + cExpData->nInputMorphogens);
      double val2 = targetSimState->getProduct(row, col, k);
      double absSub = fabs(val1 - val2) - cExpData->distErrorThreshold;
        
		  if (absSub>0){
        double d = log(1 + absSub);
        dist += d;
      }
    }
  }
    
  dist /= (cThreadData->dimX * cThreadData->dimY);

  cSharedData->reductionArray[threadIdx.x] = dist;

  __syncthreads();

  for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1) {
    if (threadIdx.x < s)
      cSharedData->reductionArray[threadIdx.x] = cSharedData->reductionArray[threadIdx.x] + 
                                                  cSharedData->reductionArray[threadIdx.x + s];
    __syncthreads();
  }

  dist = cSharedData->reductionArray[0];

  return dist;
}

  
__device__ void initializeSimState(CModelData *cModelData, CThreadData *cThreadData, 
    CSimTempData *cSimTempData, CSimState *inSimState) {
  int nInProducts = inSimState->nProducts(); // This is the number of input products (used for initialization of the state)
  int nProducts = cModelData->nProducts; // This is the number of products in the model, used during the simulation
  for (int i = cThreadData->startIndex; i < cThreadData->endIndex; ++i) {
    int row = i / cThreadData->dimX;
    int col = i % cThreadData->dimX;
    for (int k = 0; k < nInProducts; ++k) {
      double val1 = inSimState->getProduct(row, col, k);
      cSimTempData->cSimState1.setProduct(row, col, k, val1);
      cSimTempData->cSimState2.setProduct(row, col, k, val1); // Necessary because simulateExperiment only writes when conc changes
    }
    for (int k = nInProducts; k < nProducts; ++k) {
      cSimTempData->cSimState1.setProduct(row, col, k, 0.0);
      cSimTempData->cSimState2.setProduct(row, col, k, 0.0); // Necessary because simulateExperiment only writes when conc changes
    }
  }
}


__device__ double simulateExperiment(CExperimentsData *cExpData, CModelData *cModelData, 
    CSharedData *cSharedData, CThreadData *cThreadData, CSimTempData *cSimTempData, 
    CSimState *inSimState, CSimState **outSimState) {
  initializeSimState(cModelData, cThreadData, cSimTempData, inSimState);

  CSimState *simStateA = &cSimTempData->cSimState1;
  CSimState *simStateB = &cSimTempData->cSimState2;

  double change = 1.0;
  int step = 0;
  while (step < cExpData->nMaxSteps && change > cExpData->minConcChange) {
    double maxChange = 0.0;
    for (int i = cThreadData->startIndex; i < cThreadData->endIndex; ++i) {
      int row = i / cThreadData->dimX;
      int col = i % cThreadData->dimX;
      for (int k = 0; k < cModelData->nProducts; ++k)
        cSimTempData->oldConcs[blockDim.x*k + threadIdx.x] = simStateA->getProduct(row, col, k);

      computeOps(cModelData->ops, cModelData->nOps);
      computeRatios(cModelData, cThreadData, cSimTempData, row, col, simStateA);

      for (int k = 0; k < cModelData->nProducts; ++k){
        double ratio = cSimTempData->ratios[blockDim.x * k + threadIdx.x];

        if (ratio){
          double c = cSimTempData->oldConcs[blockDim.x * k + threadIdx.x] + cExpData->dt * ratio;
          if (c < cExpData->minConc)
            simStateB->setProduct(row, col, k, 0);
          else
            simStateB->setProduct(row, col, k, c);

          maxChange = fmax(maxChange, fabs(ratio));
        }
      }
    }

    cSharedData->reductionArray[threadIdx.x] = maxChange;

    __syncthreads();
      
    for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1) {
      if (threadIdx.x < s)
        cSharedData->reductionArray[threadIdx.x] = fmax(cSharedData->reductionArray[threadIdx.x], 
                                                        cSharedData->reductionArray[threadIdx.x + s]);
      __syncthreads();
    }

    change = cSharedData->reductionArray[0];

    CSimState *temp = simStateA;
    simStateA = simStateB;
    simStateB = temp;
      
    ++step;
  }

  *outSimState = simStateA; // simStateA can be simState1 or simState2

  return change;
}

  
__global__ void __launch_bounds__(NTHREADS, 1) simulateModelKernel(CExperimentsData *cExpData, 
    CModelData *cModelData, CSimTempData *cSimTempData, double *return_error) {
  CSharedData cSharedData;
  cSharedData.initialize();

  CThreadData cThreadData;
  cThreadData.initialize(cExpData);

  if (threadIdx.x == 0) {
    for (int i = 0; i < cModelData->nOps; ++i)
      cModelData->ops[i].linkFuncPointer(cSimTempData->ratios, cSimTempData->oldConcs, cSharedData.tempHill, blockDim.x);
  }

  __syncthreads();
    
  int i = 0;
  double change = 0.0;
  double error = 0.0;
  while (i < cExpData->nExperiments && error <= cModelData->maxError) {
    CSimState *inSimState = cExpData->inSimStates[i];
    CSimState *outSimState;
    change += simulateExperiment(cExpData, cModelData, &cSharedData, &cThreadData, cSimTempData, inSimState, &outSimState);

    CSimState *targetSimState = cExpData->targetSimStates[i];
    double expError = computeError(cExpData, &cSharedData, &cThreadData, outSimState, targetSimState);

    if (change > cExpData->minConcChange)
      expError += change - cExpData->minConcChange;
      
    error += expError / cExpData->nExperiments;

    ++i;
  }

  *return_error = error;
}
    
void launchKernel(CExperimentsData *cExperimentsDataDevice, CModelData *cModelDataDevice, 
    CSimTempData *cSimTempDataDevice, int nProducts, int nThreads, cudaStream_t stream, 
    double *return_error) {
  const int sharedMemorySize = sizeof(double) * (2 * NTHREADS);

  if (sharedMemorySize > MAXSHAREDMEM)
    Log::write() << "launchKernel: ERROR: Requesting " <<
      sharedMemorySize << " shared memory of " << MAXSHAREDMEM << " available. (nProducts = " << nProducts << ")" << endl;

  simulateModelKernel << < 1, NTHREADS, sharedMemorySize, stream >> > (cExperimentsDataDevice, cModelDataDevice, cSimTempDataDevice, return_error);

}

}

