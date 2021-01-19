// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#include "modelsimulator.h"

#include "simstate.h"
#include "simparams.h"
#include "simop.h"
#include "simopone.h"
#include "simopzero.h"
#include "simophillact.h"
#include "simophillrep.h"
#include "simopcopy.h"
#include "simopand.h"
#include "simopor.h"

#include "Model/model.h"
#include "Model/modelprod.h"
#include "Model/modellink.h"

namespace LoboLab {

ModelSimulator::ModelSimulator(double dift, int dx, int dy, double minConc)
  : dt_(dift), dimX_(dx), dimY_(dy), minConc_(minConc), nProducts_(0), nAllocatedProducts_(0),
    oldConcs_(NULL), ratios_(NULL),
    limits_(NULL), degradations_(NULL), nDif_(0), difProdInd_(NULL), difConsts_(NULL),
    nOps_(0), nAllocatedOps_(0), ops_(NULL) {
}

ModelSimulator::~ModelSimulator() {
  clearAll();
}

void ModelSimulator::clearAll() {
  clearLabels();
  clearProducts();
  clearOps();
}

void ModelSimulator::clearLabels() {
  labels_.clear();
  labelsInd_.clear();
  oldDiffProds_.clear();
}

void ModelSimulator::clearProducts() {
  delete [] oldConcs_;
  delete [] ratios_;
  delete [] limits_;
  delete [] degradations_;
  delete [] difProdInd_;
  delete [] difConsts_;
    
  nAllocatedProducts_ = 0;
  nDif_ = 0;
}

void ModelSimulator::clearOps() {
  deleteOps();

  delete[] ops_;

  nAllocatedOps_ = 0;
}

void ModelSimulator::deleteOps() {
  for (int i = 0; i < nOps_; ++i)
    delete ops_[i];
}

void ModelSimulator::loadModel(const Model &model, int nMorphogens) {
  clearLabels();
  QSet<int> labelSet = model.calcProductLabelsInUse(nMorphogens);
  labels_ = labelSet.toList();
  qSort(labels_);
  nProducts_ = labels_.size();

  labelsInd_ = QHash<int, int>();
  for (int i = 0; i < nProducts_; ++i)
    labelsInd_[labels_.at(i)] = i;


  if (nProducts_ > nAllocatedProducts_) {
    clearProducts();

    oldConcs_ = new double[nProducts_];
    ratios_ = new double[nProducts_];
    limits_ = new double[nProducts_];
    degradations_ = new double[nProducts_];
    difProdInd_ = new int[nProducts_];
    difConsts_ = new double[nProducts_];

    nAllocatedProducts_ = nProducts_;
  }

  QList<SimOp*> opsList; // Temporary storage for the operations

  nDif_ = 0;
  for (int i = 0; i < nProducts_; ++i) {
    // Process product constants
    ModelProd *prod = model.prodWithLabel(labels_.at(i));
    limits_[i] = prod->lim();
    degradations_[i] = prod->deg();

    // Process product diffusion
    if (prod->dif() > 0) {
      difProdInd_[nDif_] = i;
      difConsts_[nDif_] = prod->dif();
      oldDiffProds_.append(Eigen::MatrixXd());
      nDif_++;
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
    else
      opsList.append(new SimOpZero(&ratios_[i]));
  }

  int nNewOps = opsList.size();

  if (nNewOps > nAllocatedOps_) {
    clearOps();
    ops_ = new SimOp*[nNewOps];
    nAllocatedOps_ = nNewOps;
  } else
    deleteOps();

  nOps_ = nNewOps;

  for (int i = 0; i < nOps_; ++i)
    ops_[i] = opsList.at(i);
}

QList<SimOp*> ModelSimulator::createProductOps(int p, 
                              const QList<ModelLink*> &orLinks,
                              const QList<ModelLink*> &andLinks) {
  QList<SimOp*> opsList;
  bool ratiosTempUsed = false;

  // Process OR links
  int n = orLinks.size();
  for (int i = 0; i < n; ++i) {
    ModelLink *link = orLinks[i];
    if (ratiosTempUsed) {
      opsList.append(createHillOpForLink(link, &ratiosTempHill_));
      opsList.append(new SimOpOr(&ratiosTempHill_, &ratios_[p]));
    } else {
      opsList.append(createHillOpForLink(link, &ratios_[p]));
      ratiosTempUsed = true;
    }
  }

  // Process AND links
  n = andLinks.size();
  for (int i = 0; i < n; ++i) {
    ModelLink *link = andLinks[i];
    if (ratiosTempUsed) {
      opsList.append(createHillOpForLink(link, &ratiosTempHill_));
      opsList.append(new SimOpAnd(&ratiosTempHill_, &ratios_[p]));
    } else {
      opsList.append(createHillOpForLink(link, &ratios_[p]));
      ratiosTempUsed = true;
    }
  }

  if (!ratiosTempUsed) // No regulators
    opsList.append(new SimOpZero(&ratios_[p]));
    
  return opsList;
}

SimOp *ModelSimulator::createHillOpForLink(ModelLink *link, double *to) const {
  if (link->hillCoef() >= 0)
    return new SimOpHillAct(
             &oldConcs_[labelsInd_[link->regulatorProdLabel()]],
             to, link->disConst(), link->hillCoef());
  else
    return new SimOpHillRep(
             &oldConcs_[labelsInd_[link->regulatorProdLabel()]],
             to, link->disConst(), -1.0 * link->hillCoef());
}

double ModelSimulator::simulateStep(SimState *state) {
  double maxChange = 0;
  //const Eigen::MatrixXb &cells = state->cells();

  // Store old concent of diffusive products
  for (int d = 0; d < nDif_; ++d)
    oldDiffProds_[d] = state->product(difProdInd_[d]);

  for (int j=0; j < dimY_; ++j) {
    for (int i=0; i < dimX_; ++i) {
        for (int k = 0; k < nProducts_; ++k)
          oldConcs_[k] = state->product(k)(i,j);

        // Compute operations (saved in ratios)
        for (int k = 0; k < nOps_; ++k) {
          ops_[k]->compute();

        }
        // Compute ratios
        for (int k = 0; k < nProducts_; ++k) {
          ratios_[k] = limits_[k] * ratios_[k] -
            degradations_[k] * oldConcs_[k];
        }
        // Add diffusion
        for (int d = 0; d < nDif_; ++d) {
          const Eigen::MatrixXd &prod = oldDiffProds_.at(d);

          double diffusion = 0;
          double total = 0;
          if (i>0) {
            diffusion += prod(i - 1, j);
            ++total;
          }

          if (i<dimX_-1) {
            diffusion += prod(i+1,j);
            ++total;
          }

          if (j>0) {
            diffusion += prod(i, j - 1);
            ++total;
          }

          if (j<dimY_-1) {
            diffusion += prod(i, j + 1);
            ++total;
          }

          int k = difProdInd_[d];
          if (diffusion > 0 || (total > 0 && oldConcs_[k] > 0)) {
            // See nerical recipes, p. 1027 eq (20.0.6)
            // See also Schwartz2008, p. 232
            ratios_[k] += difConsts_[d] *(diffusion - total * oldConcs_[k]);
          }
       } 

        // Update product concentration 
        for (int k = 0; k < nProducts_; ++k) {
          double ratio = ratios_[k];
          if (ratio) {
            double c = oldConcs_[k] + dt_ * ratio;

            if (c < minConc_)
              state->product(k)(i,j) = 0;
            else
              state->product(k)(i, j) = c;

            maxChange = MathAlgo::max(maxChange, qAbs(ratio));
          }
        }
      }
    }

  return maxChange;
}

void ModelSimulator::blockProductDiffusion(int label) {
  int ind = labelsInd_.value(label, nProducts_);

  if (ind < nProducts_) {
    int i = 0;
    while (i < nDif_) {
      if (difProdInd_[i] == ind) {
        --nDif_;
        for (; i < nDif_; ++i) 
          difProdInd_[i] = difProdInd_[i+1];
      }
      ++i;
    }
  }
}

void ModelSimulator::blockProductProduction(int label) {
  int ind = labelsInd_.value(label, nProducts_);

  if (ind < nProducts_)
    limits_[ind] = 0;
  }

}