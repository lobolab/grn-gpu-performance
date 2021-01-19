// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#pragma once

#include <QList>

namespace LoboLab {

class SearchParams;
class SimParams;
class Search;
class Generation;
class Deme;
class Individual;
class DB;

class SearchAlgoDetCrowd {
 public:
  SearchAlgoDetCrowd(Deme *deme, Search *s);
  ~SearchAlgoDetCrowd(void);
  
  inline Deme *deme() const { return deme_; }
  inline Generation *currentGeneration() const { return generation_; }

  QList<Individual*> calcInitialPopulation();
  const QList<Individual*> &reproduce();
  void chooseNextGeneration();

 private:
  Individual *newRandIndividual(int nProducts, int nMorphogens) const;
  
  Search *search_;
  Deme *deme_;

  SearchParams *searchParams_;
  SimParams *simParams_;
  Generation *generation_;
  QList<Individual*> children_;
  int populationSize_;
  int *randPopulationInd_;
  int nMorphogens_;
  int nMinProducts_; 
  int nMaxProducts_;
  int nMaxLinks_;
};

} // namespace LoboLab
