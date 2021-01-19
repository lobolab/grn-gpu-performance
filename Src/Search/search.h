// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#pragma once

#define SearchAlgo SearchAlgoDetCrowd

#include "DB/dbelementdata.h"
#include "searchexperiment.h"

#include <QList>
#include <QDateTime>
#include <QElapsedTimer>

namespace LoboLab {

class SearchExperiment;
class Deme;
class Generation;
class Individual;
class Experiment;
class SearchParams;
class SimParams;
class ErrorCalculator;
class SearchAlgo;

class Search : public DBElement {
 public:
  Search(int id, DB *db, bool loadEvolution = true);
  ~Search();

  inline const QString &name() const { return name_; }
  inline const QDateTime &startDatetime() const { return startDatetime_; }
  inline const QDateTime &endDatetime() const { return endDatetime_; }
  inline SearchParams *searchParams() const { return searchParams_; }
  inline SimParams *simParams() const { return simParams_; }
  inline unsigned int randSeed() const { return randSeed_; }

  inline void setName(const QString &newName) { name_ = newName; }
  inline void setRandSeed(unsigned int randSeed) { randSeed_ = randSeed; }

  inline void setStartDatetime(const QDateTime &startDatetime) 
  {startDatetime_ = startDatetime;}
  inline void setEndDatetime(const QDateTime &endDatetime)
  {endDatetime_ = endDatetime;}

  inline int nExperiments() const { return searchExperiments_.size(); }
  inline Experiment *experiment(int i) const { 
    return searchExperiments_.at(i)->experiment(); 
  }
  void addExperiment(Experiment *experiment);
  bool removeExperiment(Experiment *experiment);

  inline const QList<Deme*> &demes() const {return demes_;}
  
  // Functions used during evolution by the search algorithm
  void addNewIndividual(Individual *ind);
  void removeIndividual(Individual *ind);

  void runEvolution(ErrorCalculator *errorCalculator);

  inline virtual int id() const { return ed_.id(); };
  virtual int submit(DB *db);
  int submitDemes(DB *db);
  int submitParetoFront(DB *db);
  int submitEvolution(DB *db);
  int submitShallow(DB *db);
  virtual bool erase();

 private:
  Search(const Search &source, bool maintainId = true);
  Search &operator=(const Search &source);
  void deleteAll();

  
  void calcParalEvolution(ErrorCalculator *errorCalculator);
  void processMigrationsAndReproduce(int iDeme, double meanGeneration, 
    const QList<SearchAlgo*> &searchAlgors, ErrorCalculator *errorCalculator,
    int *nMigrations, int *migrationStatus);
  void migrateIndividuals(Generation *gen1, Generation *gen2);
  void createMigrationPartners(int *migrationStatus);
  void releaseWaitingDemes(const QList<SearchAlgo*> &searchAlgors, 
    ErrorCalculator *errorCalculator, int *migrationStatus);

  void recalculateParetoFront();

  void load(bool loadEvolutionData);
  void loadSearchExperiments();
  QHash<int, Individual*> loadIndividuals();
  void loadDemes(const QHash<int, Individual*> &individualsIdMap);

  QString name_;
  SearchParams *searchParams_;
  SimParams *simParams_;
  unsigned int randSeed_;
  QDateTime startDatetime_;
  QDateTime endDatetime_;

  QList<SearchExperiment*> searchExperiments_;
  QList<Deme*> demes_;
  QList<Individual*> individuals_;
  QList<Individual*> newIndividuals_; // Temporary storage of selected 
                                      // new individuals
  QList<Individual*> paretoFront_; // Pareto front individuals 
                                   // Sorted by complexity and then by error
  QList<Individual*> oldParetoFrontInds_; // Individuals that were in the 
                                          // pareto front.
  
  DBElementData ed_;

// Persistence fields
 public:
  enum {
    FName = 1,
    FSearchParams,
    FSimParams,
    FRandSeed,
    FStartDatetime,
    FEndDatetime
  };
};

} // namespace LoboLab
