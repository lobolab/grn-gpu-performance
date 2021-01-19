// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#include "evaluatormorphogens.h"
#include "search.h"
#include "searchparams.h"
#include "Simulator/simparams.h"
#include "Simulator/simulator.h"
#include "Simulator/simstate.h"
#include "Experiment/experiment.h"
#include "Model/model.h"

#include "DB/db.h"
#include "Common/log.h"
#include "Common/mathalgo.h"

namespace LoboLab {

EvaluatorMorphogens::EvaluatorMorphogens(const Search &search) 
  : search_(search) {
  SimParams *simParams = search.simParams();
  nSims_ = simParams->nSims;
  manipulationSim_ = simParams->NumSimSteps > 0;
  minConcChange_ = simParams->minConcChange;
  errorPrecision_ = simParams->errorPrecision;

  simulator_ = new Simulator(search);

  nExperiments_ = search.nExperiments();
  for (int i = 0; i < nExperiments_; ++i) {
    Experiment *exp = search.experiment(i);
    SimState *state = new SimState(simParams->size,
                                   simParams->nInputMorphogens, simParams->nTargetMorphogens, simParams->distErrorThreshold);
    state->loadMorphologyImage(*exp->outputMorphology(), simParams->nTargetMorphogens);
    targetStates_.append(state);
  }
}

EvaluatorMorphogens::~EvaluatorMorphogens() {
  int n = targetStates_.size();
  for (int i = 0; i < n; ++i) {
    delete targetStates_.at(i);
  }

  delete simulator_;
}

double EvaluatorMorphogens::evaluate(const Model &model, double maxError) {
  double error = 0.0;
  
  simulator_->loadModel(&model);
  
  int i = 0;
  double change = 0;
  while (i < nExperiments_ && error <= maxError) {
    simulator_->loadExperiment(search_.experiment(i));
    simulator_->initialize();
    change += simulator_->simulate(search_.simParams()->NumSimSteps);
	  double expError = simulator_->simulatedState().negativeLogLikelihood(
		                          								      *targetStates_.at(i));

    if (change > minConcChange_)
      expError += (change - minConcChange_);

    error += expError / nExperiments_;

    ++i;
  }
  
  if (error > maxError)
    error = 100;
  else
    error = MathAlgo::ceilS(error, errorPrecision_);

  return error;
}

}