// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#include "cmodelsimulator.h"

#include "Model/model.h"
#include "Model/modelprod.h"
#include "Model/modellink.h"
#include "Simulator/simulatorconfig.h"
#include "Simulator/simstate.h"
#include "csimstate.cuh"
#include "csimop.cuh"
#include "cerrorhandle.h"
#include "modelsimulatordevice.cuh"
#include "cexperimentsdata.h"
#include "cmodeldata.h"
#include "csimtempdata.h"

#include <QElapsedTimer>

namespace LoboLab {

CModelSimulator::CModelSimulator(int dimX, int dimY, int nTotalMorphogens, 
                                  cudaStream_t stream)
  : dimX_(dimX), dimY_(dimY), stream_(stream), nAllocatedProducts_(0), 
    nAllocatedOps_(0) {

  cErrorHandle(cudaHostAlloc(&cModelDataHostDevice_, sizeof(CModelData), cudaHostAllocMapped));

  cErrorHandle(cudaMallocHost(&cSimTempDataHost_, sizeof(CSimTempData)));
  cSimTempDataHost_->cSimState1.initialize(dimX_, dimY_, 0);
  cSimTempDataHost_->cSimState2.initialize(dimX_, dimY_, 0);
  
  cErrorHandle(cudaMalloc(&cSimTempDataDevice_, sizeof(CSimTempData)));

  // Mapped memory to avoid a copy of a single double; improves performance notably
  cErrorHandle(cudaHostAlloc(&cErrorHostDevice_, sizeof(double), cudaHostAllocMapped));

  // Default allocation
  // Here to avoid mem allocations during first kernels
  const int initialNumProducts = 4;
  const int initialNumOps = 10;

  allocateProducts(initialNumProducts);
  allocateOps(initialNumOps);
  
}

CModelSimulator::~CModelSimulator() {
  clearAll();

  cErrorHandle(cudaFreeHost(cModelDataHostDevice_));

  cErrorHandle(cudaFree(cSimTempDataDevice_));
  cSimTempDataHost_->cSimState1.freeMatrix();
  cSimTempDataHost_->cSimState2.freeMatrix();
  cErrorHandle(cudaFreeHost(cSimTempDataHost_));

  cErrorHandle(cudaFreeHost(cErrorHostDevice_));
}

void CModelSimulator::clearAll() {
  clearLabels();
  clearProducts();
  clearOps();
}


void CModelSimulator::clearLabels() {
  labels_.clear();
  labelsInd_.clear();
}

void CModelSimulator::clearProducts() {
  cErrorHandle(cudaFreeHost(cModelDataHostDevice_->limits));
  cErrorHandle(cudaFreeHost(cModelDataHostDevice_->degradations));
  cErrorHandle(cudaFreeHost(cModelDataHostDevice_->difProdInd));
  cErrorHandle(cudaFreeHost(cModelDataHostDevice_->difConsts));

  cErrorHandle(cudaFree(cSimTempDataHost_->ratios));
  cErrorHandle(cudaFree(cSimTempDataHost_->oldConcs));
    
  nAllocatedProducts_ = 0;
}

void CModelSimulator::allocateProducts(int nProducts) {
  cErrorHandle(cudaHostAlloc(&cModelDataHostDevice_->limits, sizeof(double)*nProducts, cudaHostAllocMapped));
  cErrorHandle(cudaHostAlloc(&cModelDataHostDevice_->degradations, sizeof(double)*nProducts, cudaHostAllocMapped));
  cErrorHandle(cudaHostAlloc(&cModelDataHostDevice_->difProdInd, sizeof(double)*nProducts, cudaHostAllocMapped));
  cErrorHandle(cudaHostAlloc(&cModelDataHostDevice_->difConsts, sizeof(double)*nProducts, cudaHostAllocMapped));

  cSimTempDataHost_->cSimState1.resize(nProducts);
  cSimTempDataHost_->cSimState2.resize(nProducts);

  cErrorHandle(cudaMalloc(&cSimTempDataHost_->ratios, sizeof(double)*NTHREADS*nProducts));
  cErrorHandle(cudaMalloc(&cSimTempDataHost_->oldConcs, sizeof(double)*NTHREADS*nProducts));

  cErrorHandle(cudaMemcpyAsync(cSimTempDataDevice_, cSimTempDataHost_, sizeof(CSimTempData), cudaMemcpyHostToDevice, stream_));

  nAllocatedProducts_ = nProducts;
}


void CModelSimulator::clearOps() {
  cErrorHandle(cudaFreeHost(cModelDataHostDevice_->ops));
  nAllocatedOps_ = 0;
}

void CModelSimulator::allocateOps(int nOps) {
  cErrorHandle(cudaHostAlloc(&cModelDataHostDevice_->ops, sizeof(CSimOp) * nOps, cudaHostAllocMapped));
  nAllocatedOps_ = nOps;
}

void CModelSimulator::loadModel(const Model &model, int nMorphogens) {
  clearLabels();
  QSet<int> labelSet = model.calcProductLabelsInUse(nMorphogens);
  labels_ = labelSet.toList();
  qSort(labels_);
  int nProducts = labels_.size();
  cModelDataHostDevice_->nProducts = nProducts;

  labelsInd_ = QHash<int, int>();
  for (int i = 0; i < nProducts; ++i)
    labelsInd_[labels_.at(i)] = i;

  // Notice that reserving and freeing memory in the device prevents kernel 
  // concurrency.
  if (nProducts > nAllocatedProducts_) {
    clearProducts();
    allocateProducts(nProducts);
  }

  QList< CSimOp > opsList; // Temporary storage for the operations

  cModelDataHostDevice_->nDif = 0;
  for (int i = 0; i < nProducts; ++i) {
    // Process product constants
    ModelProd *prod = model.prodWithLabel(labels_.at(i));
    cModelDataHostDevice_->limits[i] = prod->lim();
    cModelDataHostDevice_->degradations[i] = prod->deg();

    // Process product diffusion
    if (prod->dif() > 0) {
      cModelDataHostDevice_->difProdInd[cModelDataHostDevice_->nDif] = i;
      cModelDataHostDevice_->difConsts[cModelDataHostDevice_->nDif] = prod->dif();
      cModelDataHostDevice_->nDif++;
    }

    // Process product links
    QList<ModelLink*> links = model.linksToLabel(labels_.at(i));
    QList<ModelLink*> orLinks;
    QList<ModelLink*> andLinks;

    // Categorize links
    bool anyActivator = false;
    int n = links.size();
    for (int j = 0; j < n; ++j) {
      ModelLink *link = links[j];
      if (labelSet.contains(link->regulatorProdLabel())) {
        if (link->isAndReg())
          andLinks.append(link);
        else
          orLinks.append(link);
      }
      if (link->hillCoef() >= 0)
        anyActivator = true;
    }
    if (anyActivator)
      opsList.append(createProductOps(i, orLinks, andLinks));
    else {
      opsList.append(CSimOp(CSimOp::OpZero, i));
    }
  }

  cModelDataHostDevice_->nOps = opsList.size();
  
  if (cModelDataHostDevice_->nOps > nAllocatedOps_) {
    clearOps();
    allocateOps(cModelDataHostDevice_->nOps);
  }

  for (int i = 0; i < cModelDataHostDevice_->nOps; ++i)
    cModelDataHostDevice_->ops[i] = opsList.at(i);
}

QList<CSimOp> CModelSimulator::createProductOps(int p, const QList<ModelLink*> &orLinks, const QList<ModelLink*> &andLinks) {
  QList<CSimOp> opsList;
  bool ratiosTempUsed = false;

  // Process OR links
  int n = orLinks.size();
  for (int i = 0; i < n; ++i) {
    ModelLink *link = orLinks[i];
    if (ratiosTempUsed) {
      opsList.append(createHillOpForLink(link, -1.0));
      opsList.append(CSimOp(CSimOp::OpOr, p));
    } else {
      opsList.append(createHillOpForLink(link, p));
      ratiosTempUsed = true;
    }
  }

  // Process AND links
  n = andLinks.size();
  for (int i = 0; i < n; ++i) {
    ModelLink *link = andLinks[i];
    if (ratiosTempUsed) {
      opsList.append(createHillOpForLink(link, -1.0));
      opsList.append(CSimOp(CSimOp::OpAnd, p));
    } else {
      opsList.append(createHillOpForLink(link, p));
      ratiosTempUsed = true;
    }
  }

  if (!ratiosTempUsed) // No regulators
    opsList.append(CSimOp(CSimOp::OpZero, p));
  return opsList;
}

CSimOp CModelSimulator::createHillOpForLink(ModelLink *link, double to) const {
  if (link->hillCoef() >= 0) {
    return CSimOp(CSimOp::OpHillAct, to, labelsInd_[link->regulatorProdLabel()], link->disConst(), link->hillCoef());
  } else {
    return CSimOp(CSimOp::OpHillRep, to, labelsInd_[link->regulatorProdLabel()], link->disConst(), -1.0 * link->hillCoef());
  }
}

double CModelSimulator::simulateExperiment(CExperimentsData *cExperimentsDataDevice_, int nGPUThreads_, double maxError) {
  cModelDataHostDevice_->maxError = maxError;
  //cErrorHandle(cudaMemcpyAsync(cModelDataDevice_, cModelDataHost_, sizeof(CModelData), cudaMemcpyHostToDevice, stream_));

  // default error, in case of problems
  *cErrorHostDevice_ = -1;

  launchKernel(cExperimentsDataDevice_, cModelDataHostDevice_, cSimTempDataDevice_, cModelDataHostDevice_->nProducts, nGPUThreads_, stream_, cErrorHostDevice_);

  cErrorHandle(cudaPeekAtLastError());
  cErrorHandle(cudaStreamSynchronize(stream_));

  return *cErrorHostDevice_;
}

}


