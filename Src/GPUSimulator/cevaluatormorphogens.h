// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#pragma once

#include <QList>
#include <QMutex>

#include <cuda_runtime.h>

namespace LoboLab {

class Search;
class Model;
class SimStateMap;
class SimState;
struct CSimState;
class CModelSimulator;
struct CExperimentsData;

class CEvaluatorMorphogens {
  public:
    explicit CEvaluatorMorphogens(int threadId, const Search &search);
    ~CEvaluatorMorphogens();

    double evaluate(const Model &model, double maxError);

  private:
    struct GPUExperimentDataset {
      CSimState **inStatesHost;
      CSimState **inStatesDevice;
      CSimState **targetStatesHost;
      CSimState **targetStatesDevice;
      CSimState **inStatesDeviceArray;
      CSimState **targetStatesDeviceArray;
      QList<double*> tempCpuMatrixList;
      CExperimentsData * cExperimentsDataHost;
      CExperimentsData * cExperimentsDataDevice;
    };

    void createGPUExperimentDatasets(int threadId, const Search &search);
    void createGPUExperimentDataset(const Search &search);
    void deleteGPUExperimentDatasets();

    CSimState *copySimState(const SimState *simState, GPUExperimentDataset *data);
    CSimState *copyCSimStateToDevice(CSimState *cSimStateHost);


    CModelSimulator *cModelSimulator_;

    static QVector<GPUExperimentDataset *> gpuExperimentDatasets_;
    static int nInstances_;
    static QMutex mutex_;

    int gpuDeviceId_;
    cudaStream_t stream_;

    int nMorphogens_;
    int errorPrecision_;
    int nExperiments_;
    int nGPUThreads_;
};
}
