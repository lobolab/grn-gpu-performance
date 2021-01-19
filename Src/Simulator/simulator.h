// Copyright (c) Lobo Lab (lobo@umbc.edu)
// All rights reserved.

#pragma once

#include "modelsimulator.h"
#include "simstate.h"

namespace LoboLab {

class Experiment;
class Search;

class Simulator {
 public:
  Simulator(const Search &search);
  ~Simulator();

  inline const Model *model() const { return model_; }
  inline const Experiment *experiment() const { return experiment_; }
    
  inline SimState &initialState() { return initialState_; }
  inline SimState &simulatedState() { return simulatedState_; }
  inline const SimState &initialState() const { return initialState_; }
  inline const SimState &simulatedState() const { return simulatedState_; }

  inline int nProducts() const { return modelSimulator_.nProducts(); }

  inline QSize domainSize() const { return modelSimulator_.domainSize(); }
  inline int currStep() const { return currStep_; }
  inline double currT() const { return currStep_ * modelSimulator_.dt(); }

  inline void blockProdDif(int label) { 
    modelSimulator_.blockProductDiffusion(label); 
  }
  inline void blockProdProd(int label) { 
    modelSimulator_.blockProductProduction(label); 
  }
  
  void loadModel(const Model *model);
  void loadExperiment(const Experiment *exp);
  void initialize();
  
  double simulate(int nSteps);

 private:
  const Search &search_;
  const Model *model_;
  const Experiment *experiment_;
  
  ModelSimulator modelSimulator_;
  
  SimState initialState_; 
  SimState simulatedState_;

  int currStep_;
  const double minConcChange_;
};

} // namespace LoboLab
