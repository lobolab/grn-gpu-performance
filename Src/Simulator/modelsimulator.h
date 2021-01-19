// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#pragma once

#include <QList>
#include <QHash>
#include <QSize>

#include "Common/mathalgo.h"

namespace LoboLab {

class SimState;
class SimOp;
class Model;
class ModelLink;

class ModelSimulator {
 public:
  ModelSimulator(double dt, int dimX, int dimY, double minConc);
  ~ModelSimulator();
  
  void loadModel(const Model &model, int nMorphogens);
  double simulateStep(SimState *state);

  void blockProductDiffusion(int label);
  void blockProductProduction(int label);

  inline int nProducts() const { return nProducts_; }
  inline QList<int> productLabels() const { return labels_; }

  inline double dt() const { return dt_; }
  inline QSize domainSize() const { return QSize(dimX_, dimY_); }

  void clearAll();
  void clearLabels();
  void clearProducts();
  void clearOps();
  void deleteOps();

 private:
  ModelSimulator(const ModelSimulator &source);
  ModelSimulator &operator=(const ModelSimulator &source);
  QList<SimOp*> createProductOps(int p, const QList<ModelLink*> &orLinks,
                                 const QList<ModelLink*> &andLinks);
  SimOp *createHillOpForLink(ModelLink *link, double *to) const;


  QList<int> labels_;
  QHash<int, int> labelsInd_;

  double dt_;
  int dimX_;
  int dimY_;
  double minConc_;
  int nProducts_;
  int nAllocatedProducts_;

  double *oldConcs_;
  double *ratios_;
  double ratiosTempHill_;

  double *limits_;
  double *degradations_;

  int nDif_;
  int *difProdInd_;
  double *difConsts_;
  QList<Eigen::MatrixXd> oldDiffProds_;

  int nOps_;
  int nAllocatedOps_;
  SimOp **ops_;
};

} // namespace LoboLab
