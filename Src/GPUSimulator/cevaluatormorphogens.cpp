// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#include "cevaluatormorphogens.h"

#include "Search/search.h"
#include "Search/searchparams.h"
#include "Simulator/simparams.h"

#include "Simulator/simstate.h"
#include "Experiment/experiment.h"
#include "Experiment/morphologyimage.h"
#include "Model/model.h"

#include "cerrorhandle.h"
#include "csimstate.cuh"
#include "cmodelsimulator.h"
#include "DB/db.h"

#include "Common/log.h"
#include "Common/mathalgo.h"

#include <QImage>
#include "time.h"
#include "cexperimentsdata.h"



namespace LoboLab {

QVector< CEvaluatorMorphogens::GPUExperimentDataset *> CEvaluatorMorphogens::gpuExperimentDatasets_;
int CEvaluatorMorphogens::nInstances_ = 0;
QMutex CEvaluatorMorphogens::mutex_;

CEvaluatorMorphogens::CEvaluatorMorphogens(int threadId, const Search &search) {
  clock_t t_start;
  t_start = clock();

  nExperiments_ = search.nExperiments();

  SimParams *simParams = search.simParams();
  errorPrecision_ = simParams->errorPrecision;
  nGPUThreads_ = search.simParams()->nGPUThreads;
  nMorphogens_ = search.simParams()->nInputMorphogens +
    search.simParams()->nTargetMorphogens;

  createGPUExperimentDatasets(threadId, search);

  cModelSimulator_ = new CModelSimulator(search.simParams()->size.width(),
    search.simParams()->size.height(), nMorphogens_, stream_);
    
    
  // Synchronize initialization, so it doesn't overlap with kernel launches
  // from other threads.
  cErrorHandle(cudaStreamSynchronize(stream_));
}

CEvaluatorMorphogens::~CEvaluatorMorphogens() {
  delete cModelSimulator_;
  deleteGPUExperimentDatasets();
}


void CEvaluatorMorphogens::createGPUExperimentDatasets(int threadId, const Search &search) {
  int nDevices;
  cudaGetDeviceCount(&nDevices);

  if (search.simParams()->multiGPU) {
    gpuDeviceId_ = threadId % nDevices;
    printf("GPU Thread %d: Found %d GPUs. Using GPU %d.\n", threadId, nDevices, gpuDeviceId_);
    cudaSetDevice(gpuDeviceId_);
  } else {
    gpuDeviceId_ = 0;
    printf("GPU Thread %d: Found %d GPUs. Using default GPU.\n", threadId, nDevices);
  }

  cErrorHandle(cudaStreamCreate(&stream_));

  mutex_.lock();
  nInstances_++;
    
  if (gpuExperimentDatasets_.empty()) {
    gpuExperimentDatasets_.fill(NULL, nDevices);
    createGPUExperimentDataset(search);
  } else if (!gpuExperimentDatasets_[gpuDeviceId_]) {
    createGPUExperimentDataset(search);
  }

  mutex_.unlock();
}

void CEvaluatorMorphogens::createGPUExperimentDataset(const Search &search) {
  GPUExperimentDataset * data = new GPUExperimentDataset;

  cErrorHandle(cudaMallocHost(&data->cExperimentsDataHost, sizeof(CExperimentsData)));
    
  SimParams *simParams = search.simParams();
  data->cExperimentsDataHost->nExperiments = search.nExperiments();
  data->cExperimentsDataHost->nMaxSteps = simParams->NumSimSteps;
  data->cExperimentsDataHost->minConc = simParams->minConc;
  data->cExperimentsDataHost->minConcChange = simParams->minConcChange;
  data->cExperimentsDataHost->dt = simParams->dt;
  data->cExperimentsDataHost->distErrorThreshold = simParams->distErrorThreshold;
  data->cExperimentsDataHost->nInputMorphogens = simParams->nInputMorphogens;
  data->cExperimentsDataHost->nTargetMorphogens = simParams->nTargetMorphogens;

  data->inStatesHost = new CSimState*[nExperiments_];
  cErrorHandle(cudaMallocHost(&data->inStatesDevice, sizeof(CSimState*) * nExperiments_));
  data->targetStatesHost = new CSimState*[nExperiments_];
  cErrorHandle(cudaMallocHost(&data->targetStatesDevice, sizeof(CSimState*) * nExperiments_));

  for (int i = 0; i < nExperiments_; ++i) {
    Experiment *exp = search.experiment(i);

    //*************   Load Initial State   *************************//

    SimState *state = new SimState(simParams->size, simParams->nInputMorphogens, simParams->nTargetMorphogens, simParams->distErrorThreshold);
    state->initialize(simParams->nInputMorphogens + simParams->nTargetMorphogens);
    state->loadMorphologyImage(*exp->inputMorphology(), simParams->nInputMorphogens);

    // TODO: replace states allocations with cudaMallocPitch. See https://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#programming-model
    data->inStatesHost[i] = copySimState(state, data);
    data->inStatesDevice[i] = copyCSimStateToDevice(data->inStatesHost[i]);

    delete state;

    //***************  Load Target State   ****************************// 

    SimState *tstate = new SimState(simParams->size, simParams->nInputMorphogens, simParams->nTargetMorphogens, simParams->distErrorThreshold);
    tstate->loadMorphologyImage(*exp->outputMorphology(), simParams->nTargetMorphogens);

    data->targetStatesHost[i] = copySimState(tstate, data);
    data->targetStatesDevice[i] = copyCSimStateToDevice(data->targetStatesHost[i]);

    delete tstate;
  }

  cErrorHandle(cudaMalloc(&data->inStatesDeviceArray, sizeof(CSimState*)*nExperiments_));
  cErrorHandle(cudaMemcpyAsync(data->inStatesDeviceArray, data->inStatesDevice, sizeof(CSimState*)*nExperiments_, cudaMemcpyHostToDevice, stream_));

  cErrorHandle(cudaMalloc(&data->targetStatesDeviceArray, sizeof(CSimState*)*nExperiments_));
  cErrorHandle(cudaMemcpyAsync(data->targetStatesDeviceArray, data->targetStatesDevice, sizeof(CSimState*)*nExperiments_, cudaMemcpyHostToDevice, stream_));

  data->cExperimentsDataHost->inSimStates = data->inStatesDeviceArray;
  data->cExperimentsDataHost->targetSimStates = data->targetStatesDeviceArray;

  cErrorHandle(cudaMalloc(&data->cExperimentsDataDevice, sizeof(CExperimentsData)));
  cErrorHandle(cudaMemcpyAsync(data->cExperimentsDataDevice, data->cExperimentsDataHost, sizeof(CExperimentsData), cudaMemcpyHostToDevice, stream_));

  gpuExperimentDatasets_[gpuDeviceId_] = data;
}

// Here in order to avoid CSimState to call Eigen functions, and hence avoid compile the Eigen library by nvcc.
CSimState *CEvaluatorMorphogens::copySimState(const SimState *simState, GPUExperimentDataset *data) {
  int width = simState->size().width();
  int height = simState->size().height();
  int nProducts = simState->nProducts();

  CSimState *cSimState;
  cErrorHandle(cudaMallocHost(&cSimState, sizeof(CSimState)));
  cSimState->initialize(width, height, nProducts);

  int totalSize = width * height * simState->nProducts();
  double *cpuMatrix;
  cErrorHandle(cudaMallocHost(&cpuMatrix, sizeof(double)*totalSize));
  data->tempCpuMatrixList.append(cpuMatrix);

  for (int i = 0; i < width; i++) {
    for (int j = 0; j < height; j++) {
      for (int k = 0; k < nProducts; k++) {
        int index = k * height * width + i * height + j;
        Q_ASSERT(index >= 0 && index < totalSize);
        cpuMatrix[index] = simState->product(k)(i, j);
      }
    }
  }

  cErrorHandle(cudaMemcpyAsync(cSimState->matrix_, cpuMatrix, sizeof(double)*totalSize, 
    cudaMemcpyHostToDevice, stream_));

  return cSimState;
}

CSimState *CEvaluatorMorphogens::copyCSimStateToDevice(CSimState *cSimStateHost) {
  CSimState *cSimStateDevice;
  cErrorHandle(cudaMalloc(&cSimStateDevice, sizeof(CSimState)));
  cErrorHandle(cudaMemcpyAsync(cSimStateDevice, cSimStateHost, sizeof(CSimState),
    cudaMemcpyHostToDevice, stream_));

  return cSimStateDevice;
}

void CEvaluatorMorphogens::deleteGPUExperimentDatasets() {
  mutex_.lock();
  if (--nInstances_ == 0) {
    int nDatasets = gpuExperimentDatasets_.size();
    for (int iDataSet = 0; iDataSet < nDatasets; ++iDataSet) {
      GPUExperimentDataset *data = gpuExperimentDatasets_[iDataSet];
      if (data) {
        cErrorHandle(cudaFree(data->inStatesDeviceArray));
        cErrorHandle(cudaFree(data->targetStatesDeviceArray));

        for (int i = 0; i < nExperiments_; ++i) {
          cErrorHandle(cudaFree(data->inStatesDevice[i]));
          data->inStatesHost[i]->freeMatrix();
          cErrorHandle(cudaFreeHost(data->inStatesHost[i]));
          cErrorHandle(cudaFree(data->targetStatesDevice[i]));
          data->targetStatesHost[i]->freeMatrix();
          cErrorHandle(cudaFreeHost(data->targetStatesHost[i]));
        }

        cErrorHandle(cudaFreeHost(data->inStatesDevice));
        delete[] data->inStatesHost;
        cErrorHandle(cudaFreeHost(data->targetStatesDevice));
        delete[] data->targetStatesHost;

        for (int i = 0; i < data->tempCpuMatrixList.size(); ++i)
          cErrorHandle(cudaFreeHost(data->tempCpuMatrixList[i]));

        delete data;
        gpuExperimentDatasets_[iDataSet] = NULL;
      }
    }

    gpuExperimentDatasets_.clear();
  }

  mutex_.unlock();

  cudaStreamDestroy(stream_);
}

double CEvaluatorMorphogens::evaluate(const Model &model, double maxError) {
  double error = 0.0;
  cModelSimulator_->loadModel(model, nMorphogens_);

	error = cModelSimulator_->simulateExperiment(gpuExperimentDatasets_[gpuDeviceId_]->cExperimentsDataDevice, nGPUThreads_, maxError);

  if (error > maxError)
    error = 100;
  else
    error = MathAlgo::ceilS(error, errorPrecision_);

  return error;
}

}
