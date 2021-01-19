// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#include "generationindividual.h"
#include "generation.h"
#include "Model/model.h"

namespace LoboLab {

GenerationIndividual::GenerationIndividual(Generation *g, Individual *i)
  : generation_(g), individual_(i), rank_(-1), crowdDist_(-1),
    ed_("GenerationIndividual") {
  Q_ASSERT(individual_);
}

GenerationIndividual::GenerationIndividual(Generation *g, Individual *i,
    const DBElementData &ref)
  : generation_(g), individual_(i), rank_(-1), crowdDist_(-1),
    ed_("GenerationIndividual", ref) {
  Q_ASSERT(individual_);
  load();
}

GenerationIndividual::~GenerationIndividual() {
}

GenerationIndividual::GenerationIndividual(const GenerationIndividual &source,
    Generation *g, bool maintainId)
  : generation_(g), individual_(source.individual_), rank_(source.rank_),
    crowdDist_(source.crowdDist_), ed_(source.ed_, maintainId) {
}

GenerationIndividual::GenerationIndividual(const GenerationIndividual &source,
    Generation *g, Individual *ind, bool maintainId)
  : generation_(g), individual_(ind), rank_(source.rank_),
    crowdDist_(source.crowdDist_), ed_(source.ed_, maintainId) {
}


// Persistence methods

void GenerationIndividual::load() {
  rank_ = ed_.loadValue(FRank).toInt();
  crowdDist_ = ed_.loadValue(FCrowdDist).toDouble();

  ed_.loadFinished();
}

int GenerationIndividual::submit(DB *db) {
  QPair<QString, DBElement*> refMember("Generation", generation_);

  int indId = individual_->id();
  if (!indId)
    indId = individual_->submit(db);

  QHash<QString, QVariant> values;
  values.insert("Individual", indId);
  if (rank_ > -1) values.insert("Rank", rank_);
  if (crowdDist_ > -1) values.insert("CrowdDist", crowdDist_);

  return ed_.submit(db, refMember, values);
}

bool GenerationIndividual::erase() {
  return ed_.erase();
}


bool genIndErrorLessThan(GenerationIndividual *i1, GenerationIndividual *i2) {
  return (i1->error() < i2->error());
}

bool genIndErrorComplexityLessThan(GenerationIndividual *i1,
                                   GenerationIndividual *i2) {
  return i1->error() < i2->error() ||
         (i1->error() == i2->error() &&
          i1->complexity() < i2->complexity());
}

bool genIndComplexityLessThan(GenerationIndividual *i1, GenerationIndividual *i2) {
  return (i1->complexity() < i2->complexity());
}

bool genIndCrowdDistGreaterThan(GenerationIndividual *i1, GenerationIndividual *i2) {
  return (i1->crowdDist() > i2->crowdDist());
  }

}