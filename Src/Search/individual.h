// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#pragma once

#include "DB/dbelementdata.h"
#include "Model/model.h"

namespace LoboLab {

class Evolution;
class Generation;
class GenerationIndividual;

class Individual : public DBElement {
  friend class Generation;

 public:
  Individual(Model *dm, int nMorphogens, const Individual *parent1 = NULL,
             const Individual *parent2 = NULL);
  Individual(int id, DB *db);
  explicit Individual(const DBElementData &ref);
  Individual(const Individual &source, bool maintainId = true);
  ~Individual();

  inline Model *model() const { return model_; }
  inline double error() const { return error_; }
  inline int complexity() const { return modelComplexity_; }
  inline double parentError() const { return parentError_; }

  inline void setError(double error) { error_ = error; }
  inline void setSimTime(double simTime) { simTime_ = simTime; }

  void clearGenerationIndividuals();

  bool dominates(const Individual *other) const;

  inline virtual int id() const { return ed_.id(); }
  virtual int submit(DB *db);
  virtual bool erase();

 private:
  Individual &operator=(const Individual &source);

  GenerationIndividual *addedToGeneration(Generation *generation);
  GenerationIndividual *addedToGeneration(Generation *generation,
                                          const DBElementData &ref);
  void load();

  Model *model_;
  int modelComplexity_;
  double error_;
  double simTime_;
  int parent1Id_;
  int parent2Id_;
  double parentError_;

  QList<GenerationIndividual*> generationIndividuals_;

  DBElementData ed_;

// Persistence fields
 public:
  enum {
    FModel = 1,
    FComplexity,
    FError,
    FSimTime,
    FParent1,
    FParent2
  };
};

bool indErrorLessThan(Individual *i1, Individual *i2);
bool indErrorComplexityLessThan(Individual *i1, Individual *i2);

} // namespace LoboLab
