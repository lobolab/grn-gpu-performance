// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#pragma once

#include <QList>
#include <QHash>
#include <QSize>

#include <cuda_runtime.h>

namespace LoboLab {
class SimState;
class SimOp;
class Model;
class ModelLink;
struct CExperimentsData;
struct CModelData;
struct CSimTempData;
struct CSimState;
class CSimOp;
  
class CModelSimulator {
  public:
	  CModelSimulator(int dimX, int dimY, int nTotalMorphogens, cudaStream_t stream);
    ~CModelSimulator();

    void loadModel(const Model &model, int nMorphogens);
    void clearAll();
    void clearLabels();
    void clearProducts();
    void allocateProducts(int nProducts);
    void clearOps();
    void allocateOps(int nOps);

    inline QList<int> productLabels() const { return labels_; }

	  double simulateExperiment(CExperimentsData *cExperimentsDataDevice_, int nGPUThreads_, double maxError);

  private:
	  CModelSimulator(const CModelSimulator &source);
   
    QList<CSimOp> createProductOps(int p, const QList<ModelLink*> &orLinks,
                                const QList<ModelLink*> &andLinks);
    CSimOp createHillOpForLink(ModelLink *link, double to) const;


    cudaStream_t stream_;

    CModelData *cModelDataHostDevice_;

    CSimTempData *cSimTempDataHost_;
    CSimTempData *cSimTempDataDevice_;

    double *cErrorHostDevice_;

    QList<int> labels_;
    QHash<int, int> labelsInd_;

    int dimX_;
    int dimY_;
    int nAllocatedProducts_;
    int nAllocatedOps_;
  };
}