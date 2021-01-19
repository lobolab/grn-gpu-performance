// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#pragma once

#include "DB/dbelementdata.h"
#include <QSize>


namespace LoboLab {

class SimParams : public DBElement {
 public:
  SimParams(int id, DB *db);
  ~SimParams();

  SimParams(const SimParams &source, bool maintainId = true);
  SimParams &operator=(const SimParams &source);
  
  inline virtual int id() const { return ed.id(); }
  virtual int submit(DB *db);
  virtual bool erase();


  QString name;
  QSize size;
  double dt;
  double dx;
  int NumSimSteps;
  int nSims;
  int nInputMorphogens;
  int nTargetMorphogens;
  int nMaxProducts;
  int nMaxLinks;
  double minConc;
  double minConcChange;
  double distErrorThreshold;
  int errorPrecision;
  int nCPUSlaves;
  int nGPUSlaves;
  int nGPUThreads;
  int nTestIndividuals;
  bool multiGPU;

 private:
  void copy(const SimParams &source);
  void load();

  DBElementData ed;

// Persistence fields
 public:
   enum {
     FName = 1,
     FSizeX,
     FSizeY,
     FDt,
     FDx,
     FNumSimSteps,
     FNumSims,
     FNumInputMorphogens,
     FNumTargetMorphogens,
     FNumMaxProducts,
     FNumMaxLinks,
     FMinConc,
     FMinConcChange,
     FDistErrorThreshold,
     FErrorPrecision,
     FNumCPUSlaves,
     FNumGPUSlaves,
     FNumGPUThreads,
     FNumTestIndividuals,
     FMultiGPU
  };
};

} // namespace LoboLab
