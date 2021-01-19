// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#include "individual.h"
#include "generation.h"
#include "search.h"
#include "Model/model.h"
#include "DB/db.h"
#include "Common/mathalgo.h"

namespace LoboLab {

Individual::Individual(Model *dm, int nMorphogens, const Individual *parent1,
                       const Individual *parent2)
  : model_(dm), 
    error_(-1),
    simTime_(-1), 
    ed_("Individual") {
  if (parent1) { 
    parent1Id_ = parent1->id();
    parentError_ = parent1->error();
  } else {
    parent1Id_ = -1;
    parentError_ = 1000;
  }

  if (parent2)
    parent2Id_ = parent2->id();
  else
    parent2Id_ = -1;

  
  modelComplexity_ = model_->calcComplexityInUse(nMorphogens);
}

Individual::Individual(const DBElementData &ref)
  : ed_("Individual", ref) {
  load();
}

Individual::Individual(int id, DB *db)
  : ed_("Individual", id, db) {
  load();
}


Individual::~Individual() {
  clearGenerationIndividuals();
  delete model_;
}

Individual::Individual(const Individual &source, bool maintainId)
  : modelComplexity_(source.modelComplexity_),
    error_(source.error_), 
    simTime_(source.simTime_), 
    ed_(source.ed_, maintainId) {
  model_ = new Model(*source.model_);

  if (maintainId) {
    parent1Id_ = source.parent1Id_;
    parent2Id_ = source.parent2Id_;
    parentError_ = source.parentError_;
  } else {
    parent1Id_ = source.id(); 
    parent2Id_ = 0;
    parentError_ = source.error();
  }
}

void Individual::clearGenerationIndividuals() {
  int n = generationIndividuals_.size();
  for (int i = 0; i < n; ++i)
    delete generationIndividuals_.at(i);

  generationIndividuals_.clear();
}

GenerationIndividual *Individual::addedToGeneration(Generation *generation) {
  GenerationIndividual *gi = new GenerationIndividual(generation, this);
  generationIndividuals_.append(gi);

  return gi;
}

GenerationIndividual *Individual::addedToGeneration(Generation *generation,
                                                    const DBElementData &ref) {
  GenerationIndividual *gi = new GenerationIndividual(generation, this, ref);
  generationIndividuals_.append(gi);

  return gi;
}

bool Individual::dominates(const Individual *other) const {
    return (error_ < other->error_ && 
            modelComplexity_ <= other->modelComplexity_) ||
           (modelComplexity_ < other->modelComplexity_ && 
            error_ <= other->error_);
}

// Persistence methods

void Individual::load() {
  model_ = new Model();
  QString modelStr = ed_.loadValue(FModel).toString();
  model_->loadFromString(modelStr);
  modelComplexity_ = model_->calcComplexity();

  error_ = ed_.loadValue(FError).toDouble();
  simTime_ = ed_.loadValue(FSimTime).toDouble();
  parent1Id_ = ed_.loadValue(FParent1).toInt();
  parent2Id_ = ed_.loadValue(FParent2).toInt();
  parentError_ = -1;

  ed_.loadFinished();
}

int Individual::submit(DB *db) {
  QHash<QString, QVariant> values;
  values.insert("Model", model_->toString());
  values.insert("Complexity", modelComplexity_);
  values.insert("Error", error_);
  values.insert("SimTime", simTime_);
  values.insert("Parent1", parent1Id_ > -1 ? parent1Id_ : QVariant());
  values.insert("Parent2", parent2Id_ > -1 ? parent2Id_ : QVariant());

  return ed_.submit(db, values, generationIndividuals_);
}

bool Individual::erase() {
  QList<DBElement*> members;

  return ed_.erase(members, generationIndividuals_);
}

bool indErrorLessThan(Individual *i1, Individual *i2) {
  return (i1->error() < i2->error());
}

bool indErrorComplexityLessThan(Individual *i1, Individual *i2) {
  return i1->error() < i2->error() ||
         (i1->error() == i2->error() &&
          i1->complexity() < i2->complexity());
}

}
