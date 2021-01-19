// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#include "searchalgodetcrowd.h"
#include "DB/db.h"
#include "search.h"
#include "searchparams.h"
#include "Simulator/simparams.h"
#include "errorcalculator.h"
#include "deme.h"
#include "generation.h"
#include "individual.h"
#include "Model/model.h"
#include "Common/log.h"
#include "Common/mathalgo.h"
#include <QDateTime>
#include <QDebug>
#include <iostream>

namespace LoboLab {

SearchAlgoDetCrowd::SearchAlgoDetCrowd(Deme *deme, Search *s)
    : search_(s), 
      deme_(deme),
      generation_(NULL) {
  searchParams_ = search_->searchParams();
  simParams_ = search_->simParams();

  nMorphogens_ = simParams_->nInputMorphogens + 
                 simParams_->nTargetMorphogens;
  nMinProducts_ = nMorphogens_ ;
  nMaxProducts_ = simParams_->nMaxProducts;
  nMaxLinks_ = simParams_->nMaxLinks;

  populationSize_ = searchParams_->demesSize;
  if (populationSize_ % 2 != 0)
    populationSize_++;

  randPopulationInd_ = new int[populationSize_];
  for (int i = 0; i < populationSize_; ++i)
    randPopulationInd_[i] = i;
}

SearchAlgoDetCrowd::~SearchAlgoDetCrowd(void) {
  delete [] randPopulationInd_;

  if (generation_) {
    int n = generation_->nIndividuals();
    for (int i = 0; i < n ; ++i)
      search_->removeIndividual(generation_->individual(i));
  }
}

QList<Individual*> SearchAlgoDetCrowd::calcInitialPopulation() {
  generation_ = deme_->createNextGeneration();

  for (int i = 0; i < populationSize_; ++i)
    generation_->addIndividual(newRandIndividual(nMinProducts_, nMorphogens_));

  return generation_->individuals();
}

const QList<Individual*> &SearchAlgoDetCrowd::reproduce() {
  // Randomize parents
  MathAlgo::shuffle(populationSize_, randPopulationInd_);

  // Create children
  for (int i = 0; i < populationSize_; i = i + 2) {
    Individual *parent1 = generation_->individual(randPopulationInd_[i]);
    Individual *parent2 = generation_->individual(randPopulationInd_[i+1]);
    Individual *child1, *child2;
    Model *childModel1, *childModel2;

    if (MathAlgo::rand100() < 75){ // Crossover
      Model::cross(parent1->model(), parent2->model(), childModel1, childModel2, 
        nMorphogens_, nMinProducts_, nMaxProducts_, nMaxLinks_);

      // Here because Individual caches the model complexity
      childModel1->mutate(nMorphogens_, nMinProducts_, nMaxProducts_, nMaxLinks_);
      childModel2->mutate(nMorphogens_, nMinProducts_, nMaxProducts_, nMaxLinks_);

      child1 = new Individual(childModel1, nMorphogens_, parent1, parent2);
      child2 = new Individual(childModel2, nMorphogens_, parent2, parent1);
    } else { // No crossover
      childModel1 = new Model(*parent1->model());
      childModel2 = new Model(*parent2->model());

      // Here because Individual caches the model complexity
      childModel1->mutate(nMorphogens_, nMinProducts_, nMaxProducts_, nMaxLinks_);
      childModel2->mutate(nMorphogens_, nMinProducts_, nMaxProducts_, nMaxLinks_);

      child1 = new Individual(childModel1, nMorphogens_, parent1);
      child2 = new Individual(childModel2, nMorphogens_, parent2);
    }

    children_.append(child1);
    children_.append(child2);
  }

  return children_;
}

void SearchAlgoDetCrowd::chooseNextGeneration() {
  // All individuals are chosen in first generation
  if(children_.isEmpty()) {
    int numInd = generation_->nIndividuals();
    for (int i=0; i < numInd; ++i)
      search_->addNewIndividual(generation_->individual(i));
  } else {
    Generation *nextGeneration = deme_->createNextGeneration();
    
    // Select new population
    for (int i = 0; i < populationSize_; i = i + 2) {
      Individual *parent1 = generation_->individual(randPopulationInd_[i]);
      Individual *parent2 = generation_->individual(randPopulationInd_[i+1]);
      Individual *child1 = children_.at(i);
      Individual *child2 = children_.at(i+1);

      if (child1->error() <= parent1->error()) { 
        nextGeneration->addIndividual(child1);
        search_->addNewIndividual(child1);
        search_->removeIndividual(parent1);
      } else {
        nextGeneration->addIndividual(parent1);
        delete child1;
      }

      if (child2->error() <= parent2->error()) { 
        nextGeneration->addIndividual(child2);
        search_->addNewIndividual(child2);
        search_->removeIndividual(parent2);
      } else {
        nextGeneration->addIndividual(parent2);
        delete child2;
      }
    }

    children_.clear();
    generation_->clearIndividuals(); // Save some memory
    generation_ = nextGeneration;
  }

  // Minimum 1 second for better log show
  generation_->setTime(1 + search_->startDatetime().secsTo(
                                              QDateTime::currentDateTimeUtc()));
  generation_->calcPopulSta();
}

Individual *SearchAlgoDetCrowd::newRandIndividual(int nProducts, int nMorphogens) const {
  Model *model = Model::createRandom(nProducts, nMorphogens);
  Individual *newInd = new Individual(model, nMorphogens);
  return newInd;
}

}