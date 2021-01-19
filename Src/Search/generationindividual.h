// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#pragma once

#include "DB/dbelementdata.h"
#include "individual.h"

namespace LoboLab {

class Generation;

class GenerationIndividual : public DBElement {
  friend class Individual;

 public:
  inline Generation *generation() const { return generation_; }
  inline Individual *individual() const { return individual_; }

  void removeFromGeneration();

  inline double error() const { return individual_->error(); }
  inline int complexity() const { return individual_->complexity(); }
  inline int rank() const { return rank_; }
  inline double crowdDist() const { return crowdDist_; }

  inline void setRank(int rank) { rank_ = rank; }
  inline void setCrowdDist(double crowdDist) { crowdDist_ = crowdDist; }

 protected:
  inline virtual int id() const { return ed_.id(); };
  virtual int submit(DB *db);
  virtual bool erase();

 private:
  GenerationIndividual(Generation *g, Individual *i);
  GenerationIndividual(Generation *g, Individual *i,
                       const DBElementData &ref);
  ~GenerationIndividual();

  GenerationIndividual(const GenerationIndividual &source,
                       Generation *g, bool maintainId = false);
  GenerationIndividual(const GenerationIndividual &source,
                       Generation *g, Individual *ind, bool maintainId = true);
  GenerationIndividual &operator=(const GenerationIndividual &source);

  void setIndividual(Individual *newIndividual);

  void load();

  Generation *generation_;
  Individual *individual_;
  int rank_;
  double crowdDist_;

  DBElementData ed_;

// Persistence fields
 public:
  enum {
    FGeneration = 1,
    FIndividual,
    FRank,
    FCrowdDist
  };
};

bool genIndErrorLessThan(GenerationIndividual *i1, GenerationIndividual *i2);
bool genIndErrorComplexityLessThan(GenerationIndividual *i1, 
                                   GenerationIndividual *i2);
bool genIndComplexityLessThan(GenerationIndividual *i1, 
                              GenerationIndividual *i2);
bool genIndCrowdDistGreaterThan(GenerationIndividual *i1, 
                                GenerationIndividual *i2);

} // namespace LoboLab
