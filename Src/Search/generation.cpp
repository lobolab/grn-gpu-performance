// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#include "generation.h"
#include "deme.h"
#include "individual.h"
#include "Model/model.h"
#include "Common/log.h"
#include <iostream>

namespace LoboLab {

Generation::Generation(Deme *d, int i)
  : deme_(d), 
    ind_(i), 
    time_(-1),   
    minFit_(1e100),
    meanFit_(0),
    maxFit_(-1),
    minComp_(1e100),
    meanComp_(0),
    maxComp_(-1),
    bestComp_(1e100),
    ed_("Generation") {
}

Generation::Generation(Deme *d,
                       const QHash<int, Individual*> &individualsIdMap, 
                       const DBElementData &ref)
    : deme_(d), 
      ed_("Generation", ref) {
  load(individualsIdMap);
}

Generation::~Generation() {}

void Generation::addIndividual(Individual *ind) {
  individuals_.append(ind);
  ind->addedToGeneration(this);
}

void Generation::calcPopulSta() {

  int n = individuals_.size();
  for (int i = 0; i < n; ++i) {
    double fit = individuals_.at(i)->error();
    int comp = individuals_.at(i)->complexity();

    meanFit_ += fit;
    if (fit < minFit_) {
      minFit_ = fit;
      bestComp_ = comp;
    } else if (fit == minFit_ && comp < bestComp_) {
      bestComp_ = comp;
    }
    if (fit > maxFit_)
      maxFit_ = fit;

    meanComp_ += comp;
    if (comp < minComp_)
      minComp_ = comp;
    if (comp > maxComp_)
      maxComp_ = comp;
  }

  meanFit_ /= n;
  meanComp_ /= n;
  
  Log::write() << "Generation::calcPopulSta: MaxFit=" << maxFit_ <<
               " MeanFit=" << meanFit_ << " MinFit=" << minFit_ << " BestComp="
               << bestComp_ << endl;
}

// Persistence methods

void Generation::load(const QHash<int, Individual*> &individualsIdMap) {
  ind_ = ed_.loadValue(FInd).toInt();
  time_ = ed_.loadValue(FTime).toInt();

  loadGenerationIndividuals(individualsIdMap);

  ed_.loadFinished();
}

void Generation::loadGenerationIndividuals(
  const QHash<int, Individual*> &individualsIdMap) {
  ed_.loadReferences("GenerationIndividual");
  while (ed_.nextReference()) {
    Individual *ind = individualsIdMap.value(
                        ed_.loadRefValue(GenerationIndividual::FIndividual).toInt());
    individuals_.append(ind);
    ind->addedToGeneration(this, ed_);
  }
}

int Generation::submit(DB *db) {
  QPair<QString, DBElement*> refMember("deme", deme_);
  
  QHash<QString, QVariant> values;
  values.insert("Ind", ind_);
  values.insert("Time", time_);
  values.insert("MinError", minFit_);
  values.insert("MeanError", meanFit_);
  values.insert("MaxError", maxFit_);
  values.insert("MinComp", minComp_);
  values.insert("MeanComp", meanComp_);
  values.insert("MaxComp", maxComp_);
  values.insert("BestComp", bestComp_);

  QHash<QString, DBElement*> members;

  return ed_.submit(db, refMember, members, values);
}

int Generation::submitWithIndividuals(DB *db) {
  QPair<QString, DBElement*> refMember("deme", deme_);
  
  QHash<QString, QVariant> values;
  values.insert("Ind", ind_);
  values.insert("Time", time_);
  values.insert("MinError", minFit_);
  values.insert("MeanError", meanFit_);
  values.insert("MaxError", maxFit_);
  values.insert("MinComp", minComp_);
  values.insert("MeanComp", meanComp_);
  values.insert("MaxComp", maxComp_);
  values.insert("BestComp", bestComp_);

  QHash<QString, DBElement*> members;

  return ed_.submit(db, refMember, members, values, individuals_);
}

bool Generation::erase() {
  return ed_.erase();
}

}