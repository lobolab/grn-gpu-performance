// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#pragma once

#include "DB/dbelementdata.h"
#include "generationindividual.h"

namespace LoboLab {

class Individual;
class Deme;

class Generation : public DBElement {
  friend class Deme;
  friend class Search;
 public:
  ~Generation();

  inline Deme *deme() const { return deme_; }

  inline int nIndividuals() const { return individuals_.size(); }
  inline Individual *individual(int i) const { return individuals_.at(i); }
  inline const QList<Individual*> &individuals() const { return individuals_; }
  inline void clearIndividuals() { individuals_.clear(); }

  inline int ind() const { return ind_; }
  inline int time() const { return time_; }
  inline void setTime(int time) { time_ = time; }

  void calcPopulSta();

  void addIndividual(Individual *ind);

  inline virtual int id() const { return ed_.id(); }
  virtual int submit(DB *db);
  virtual int submitWithIndividuals(DB *db);
  virtual bool erase();

 private:
  Generation(Deme *d, int ind);
  Generation(Deme *d, const QHash<int, Individual*> &individualsIdMap,
             const DBElementData &ref);
  Generation(const Generation &source);
  Generation &operator=(const Generation &source);

  void load(const QHash<int, Individual*> &individualsIdMap);
  void loadGenerationIndividuals(
    const QHash<int, Individual*> &individualsIdMap);

  Deme *deme_;
  int ind_;
  int time_;
  double minFit_, meanFit_, maxFit_;
  double minComp_, meanComp_, maxComp_, bestComp_;  
  QList<Individual*> individuals_;

  DBElementData ed_;

// Persistence fields
 public:
  enum {
    FDeme = 1,
    FInd,
    FTime,
    FMinError,
    FMeanError,
    FMaxError,
    FMinComp,
    FMeanComp,
    FMaxComp,
    FBestComp
  };
};

} // namespace LoboLab
