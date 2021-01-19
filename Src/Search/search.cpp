// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#include "search.h"

#include "individual.h"
#include "searchparams.h"
#include "Simulator/simparams.h"
#include "Experiment/experiment.h"
#include "deme.h"
#include "generation.h"
#include "errorcalculator.h"
#include "DB/db.h"
#include "Common/log.h"
#include "Common/mathalgo.h"


#include "searchalgodetcrowd.h"

namespace LoboLab {
  
Search::Search(int id, DB *db, bool loadEvolution)
  : ed_("Search", id, db) {
  load(loadEvolution);
}

Search::~Search() {
  deleteAll();
}

void Search::deleteAll() {
  int n = searchExperiments_.size();
  for (int i = 0; i < n; ++i)
    delete searchExperiments_.at(i);

  searchExperiments_.clear();
  
  n = paretoFront_.size();
  for (int i = 0; i < n; ++i)
    delete paretoFront_.at(i);

  n = oldParetoFrontInds_.size();
  for (int i = 0; i < n; ++i)
    delete oldParetoFrontInds_.at(i);

  n = individuals_.size();
  for (int i = 0; i < n; ++i) {
    Individual *ind = individuals_[i];
    if (!paretoFront_.contains(ind) && !oldParetoFrontInds_.contains(ind))
      delete ind;
  }
  
  paretoFront_.clear();
  oldParetoFrontInds_.clear();
  individuals_.clear();

  n = demes_.size();
  for (int i = 0; i < n; ++i)
    delete demes_.at(i);

  demes_.clear();

  delete searchParams_;
  searchParams_ = NULL;

  delete simParams_;
  simParams_ = NULL;
}

void Search::addExperiment(Experiment *experiment) {
  SearchExperiment *se = new SearchExperiment(this, experiment);
  searchExperiments_.append(se);
}

bool Search::removeExperiment(Experiment *experiment) {
  bool taken = false;
  int n = searchExperiments_.size();
  int i = 0;
  while (!taken && i < n) {
    SearchExperiment *se = searchExperiments_.at(i);
    if (se->experiment() == experiment) {
      searchExperiments_.removeAt(i);
      ed_.removeReference(se);
      taken = true;
    }

    ++i;
  }

  return taken;
}

void Search::runEvolution(ErrorCalculator *errorCalculator) {
  name_ += " - " + QDateTime::currentDateTimeUtc().toString(DB::DatetimeFormat);
  
  startDatetime_ = QDateTime::currentDateTimeUtc();
  calcParalEvolution(errorCalculator);
  endDatetime_ = QDateTime::currentDateTimeUtc();
}

void Search::calcParalEvolution(ErrorCalculator *errorCalculator) {
  // One searchAlgorithm per deme
  QElapsedTimer timer;
  QList<SearchAlgo*> searchAlgors;

  searchAlgors.reserve(searchParams_->nDemes);
  // Create initial populations
  for (int i = 0; i < searchParams_->nDemes; ++i) {
    Deme *deme = new Deme(this);
    demes_.append(deme);
    SearchAlgo *algo = new SearchAlgo(deme, this);
    searchAlgors.append(algo);
    errorCalculator->process(i, algo->calcInitialPopulation());
  }

  double meanGeneration = 0;  
  int nMigrations = 0;
  // Migration partners (-1 if no migration is neccesary)
  int *migrationStatus = new int [searchParams_->nDemes];
  for (int i = 0; i < searchParams_->nDemes; ++i)
    migrationStatus[i] = -2;
  
  // Save during evolution
  if (searchParams_->saveIndividuals > 1)
    submitDemes(ed_.db());

  // Main loop
  int lastGeneration = searchParams_->nGenerations;
  int bestPerfectComp = 1e5;
  int nDemesLeft = searchParams_->nDemes;
  timer.start();
  while (nDemesLeft > 0) {
    Log::write() << "Search::calcParalEvolution: waiting for " << nDemesLeft << 
      " demes" << endl; 
    int iDeme = errorCalculator->waitForAnyDeme();
    Log::write() << "Search::calcParalEvolution: received deme " << iDeme << 
      endl; 

    SearchAlgo *algo = searchAlgors[iDeme];
    algo->chooseNextGeneration();
    recalculateParetoFront();

    // Save pareto front
    if (searchParams_->saveIndividuals == 2) {
      algo->currentGeneration()->submit(ed_.db());
      submitParetoFront(ed_.db());
    }
    // Save all
    else if (searchParams_->saveIndividuals == 3)
      algo->currentGeneration()->submitWithIndividuals(ed_.db());
    
    int iGen = algo->currentGeneration()->ind();
    Log::write() << "Search::calcParalEvolution: " << paretoFront_.size() << 
      " (" << oldParetoFrontInds_.size() << ") Paretos after deme " << iDeme << 
      " gen " << iGen << " meanGen " << meanGeneration << endl;

    if (meanGeneration < lastGeneration) {
      // Check if we can end prematurely
      if (paretoFront_.first()->error() == 0 &&
          paretoFront_.first()->complexity() < bestPerfectComp) {
        bestPerfectComp = paretoFront_.first()->complexity();
        lastGeneration = MathAlgo::min(250 + (int)meanGeneration, 
                                       searchParams_->nGenerations);
        Log::write() << "Search::calcParalEvolution: New bestPerfectComp found " 
          "of " << bestPerfectComp << ". Last generation set to " << 
          lastGeneration << endl;
      }

      processMigrationsAndReproduce(iDeme, meanGeneration, searchAlgors, 
        errorCalculator, &nMigrations, migrationStatus);
      meanGeneration += 1.0/searchParams_->nDemes;
    } else { // End evolution
      releaseWaitingDemes(searchAlgors, errorCalculator, migrationStatus);
      delete algo;
      --nDemesLeft;
    }

	  int s = timer.elapsed() / 1000;
    int m = (s % (60 * 60)) / 60;
	  if (m >= 1) {
      Log::write() << "Search::calcParalEvolution: saving to database..." << endl;
	    submitEvolution(ed_.db());
      Log::write() << "Search::calcParalEvolution: saved to database." << endl;

      // Remove saved old individuals
      int n = oldParetoFrontInds_.size();
      for (int i = n - 1; i >= 0; --i) {
        Individual *individual = oldParetoFrontInds_.at(i);
        if (!individuals_.contains(individual)) {
          delete individual;
          oldParetoFrontInds_.removeAt(i);
        }
      }

      // Remove old individualGenerations
      n = paretoFront_.size();
      for (int i = 0; i < n; ++i)
        paretoFront_[i]->clearGenerationIndividuals();

      n = oldParetoFrontInds_.size();
      for (int i = 0; i < n; ++i)
        oldParetoFrontInds_[i]->clearGenerationIndividuals();

      // Remove old generations
      for (int i = 0; i < searchParams_->nDemes; ++i) {
        Deme *deme = demes_[i];
        while (deme->generations().size() > 1)
          delete deme->generations().takeFirst();
      }

      timer.start();
    }
  }

  delete [] migrationStatus;
}

// Process migration status
// migrationStatus_ stores the partner to migrate with or -1 if it is waiting
// for the partner or -2 if no migration is necessary.
void Search::processMigrationsAndReproduce(int iDeme, double meanGeneration, 
    const QList<SearchAlgo*> &searchAlgors, ErrorCalculator *errorCalculator,
    int *nMigrations, int *migrationStatus) {
  // Check if migration is necessary
  int iDemeStatus = migrationStatus[iDeme];

  if (iDemeStatus > -1) { // It needs to migrate
    int iDemePartner = iDemeStatus;

    // Check if the partner is ready
    if (migrationStatus[iDemePartner] == -1) { // partner waiting
      Log::write() << "Search::processMigrationsAndReproduce: Migration "
        "between demes " << iDeme << " and " << iDemePartner << endl;
      migrationStatus[iDeme] = -2;
      migrationStatus[iDemePartner] = -2;
      migrateIndividuals(searchAlgors[iDeme]->currentGeneration(), 
                         searchAlgors[iDemePartner]->currentGeneration());
      errorCalculator->process(iDemePartner, 
                                searchAlgors[iDemePartner]->reproduce());
      errorCalculator->process(iDeme, 
                                searchAlgors[iDeme]->reproduce());
    } else {
      migrationStatus[iDeme] = -1; // waiting for partner
    }

  } else if (meanGeneration >= 
             (1 + *nMigrations) * searchParams_->migrationPeriod) { 
    // New migration
    (*nMigrations)++;
    Log::write() << "Search::processMigrationsAndReproduce: Migration "
      "number " << *nMigrations << " scheduled" << endl;
    
    // Just in case, release any deme waiting from previous migration
    releaseWaitingDemes(searchAlgors, errorCalculator, migrationStatus);
    createMigrationPartners(migrationStatus);
    migrationStatus[iDeme] = -1;

  } else {// no migration necessary
      errorCalculator->process(iDeme, searchAlgors[iDeme]->reproduce());
  }
}

void Search::migrateIndividuals(Generation *gen1, Generation *gen2) {
  int demesSize = searchParams_->demesSize;
  
  int *randIds = new int[2*demesSize];
  for (int i = 0; i < 2*demesSize; ++i)
    randIds[i] = i;

  MathAlgo::shuffle(2*demesSize, randIds);
  
  QList<Individual*> individuals;
  individuals.reserve(2*demesSize);

  for (int i = 0; i < demesSize; i++)
    individuals.append(gen1->individual(i));
  
  for (int i = 0; i < demesSize; i++)
    individuals.append(gen2->individual(i));
    
  gen1->clearIndividuals();
  gen2->clearIndividuals();

  for (int i = 0; i < demesSize; ++i)
    gen1->addIndividual(individuals[randIds[i]]);

  for (int i = demesSize; i < 2*demesSize; ++i)
    gen2->addIndividual(individuals[randIds[i]]);

  delete [] randIds;
}

void Search::createMigrationPartners(int *migrationStatus) {

  // Shuffle demes
  int *randIds = new int[searchParams_->nDemes];
  for (int i = 0; i < searchParams_->nDemes; ++i)
    randIds[i] = i;
  
  MathAlgo::shuffle(searchParams_->nDemes, randIds);

  // Select random partners
  for (int i = 0; i < searchParams_->nDemes; i = i + 2) {
    int id1 = randIds[i];
    int id2 = randIds[i+1];
    migrationStatus[id1] = id2;
    migrationStatus[id2] = id1;
  }

  delete [] randIds;
}

void Search::releaseWaitingDemes(const QList<SearchAlgo*> &searchAlgors, 
                                 ErrorCalculator *errorCalculator,
                                 int *migrationStatus) {
  for (int i=0; i < searchParams_->nDemes; ++i) {
    if (migrationStatus[i] > -2) {
      if (migrationStatus[i] == -1) {
        Log::write() << "Search::releaseWaitingDemes: Releasing waiting deme " 
          << i << endl;
        errorCalculator->process(i, searchAlgors[i]->reproduce());
      } else {
        Log::write() << "Search::releaseWaitingDemes: Canceling scheduled "
          "deme " << i << endl;
      }

      migrationStatus[i] = -2;
    }
  }
}

void Search::addNewIndividual(Individual *ind) {
  individuals_.append(ind);
  newIndividuals_.append(ind);
}

void Search::removeIndividual(Individual *ind) {
  if (paretoFront_.contains(ind)) {
    paretoFront_.removeOne(ind);
    oldParetoFrontInds_.append(ind);
  }
  else if (!oldParetoFrontInds_.contains(ind))
    delete ind;
    
  individuals_.removeOne(ind);
}

// Also includes duplicates in the front
void Search::recalculateParetoFront() {
  paretoFront_.append(newIndividuals_);

  if (!paretoFront_.isEmpty()) {
    qSort(paretoFront_.begin(), paretoFront_.end(),	indErrorComplexityLessThan);

    double lastError = paretoFront_.first()->error();
    double lastComplex = paretoFront_.first()->complexity();
    int i = 1;
    while (i < paretoFront_.size()) {
      Individual *ind = paretoFront_.at(i);
      if (ind->complexity() < lastComplex) {
        lastError = ind->error();
        lastComplex = ind->complexity();
        ++i;
      } else if (ind->complexity() == lastComplex &&
                 ind->error() == lastError) { // Duplicates included
        ++i;
      } else { // The individual is not pareto
        if (newIndividuals_.contains(ind)) // The individual was never pareto
          paretoFront_.removeAt(i);
        else  // The individual was pareto
          oldParetoFrontInds_.append(paretoFront_.takeAt(i));
      }
    }
  }

  newIndividuals_.clear();
}


// Persistence methods

void Search::load(bool loadEvolution) {
  name_ = ed_.loadValue(FName).toString();
  searchParams_ = new SearchParams(ed_.loadValue(FSearchParams).toInt(),
                                  ed_.db());
  simParams_ = new SimParams(ed_.loadValue(FSimParams).toInt(),
                            ed_.db());
  startDatetime_ = QDateTime::fromString(
                    ed_.loadValue(FStartDatetime).toString(), DB::DatetimeFormat);
  endDatetime_ = QDateTime::fromString(
                  ed_.loadValue(FEndDatetime).toString(), DB::DatetimeFormat);
  randSeed_ = ed_.loadValue(FRandSeed).toUInt();

  loadSearchExperiments();
  if (loadEvolution) {
    QHash<int, Individual*> individualsIdMap = loadIndividuals();
    loadDemes(individualsIdMap);
  }

  ed_.loadFinished();
}

void Search::loadSearchExperiments() {
  ed_.loadReferences("SearchExperiment");
  while (ed_.nextReference()) {
    SearchExperiment *ele = new SearchExperiment(this, ed_);
    searchExperiments_.append(ele);
  }
}

QHash<int, Individual*> Search::loadIndividuals() {
  QHash<int, Individual*> individualsIdMap;

  ed_.loadReferencesIndirect("Individual", "GenerationIndividual",
                            "Generation", "Deme");
  while (ed_.nextReference()) {
    Individual *ele = new Individual(ed_);
    individualsIdMap.insert(ele->id(), ele);
    individuals_.append(ele);
  }

  paretoFront_ = individuals_;
  recalculateParetoFront();

  return individualsIdMap;
}

void Search::loadDemes(const QHash<int, Individual*> &individualsIdMap) {
  ed_.loadReferences("Deme");
  while (ed_.nextReference()) {
    Deme *d = new Deme(this, individualsIdMap, ed_);
    demes_.append(d);
  }
}

int Search::submit(DB *db) {
  QPair<QString, DBElement*> refMember;

  QHash<QString, QVariant> values;
  values.insert("Name", name_);
  values.insert("RandSeed", randSeed_);
  values.insert("StartDatetime", startDatetime_.toString(DB::DatetimeFormat));
  values.insert("EndDatetime", endDatetime_.toString(DB::DatetimeFormat));

  QHash<QString, DBElement*> members;
  members.insert("SearchParams", searchParams_);
  members.insert("SimParams", simParams_);

  return ed_.submit(db, refMember, members, values, searchExperiments_, demes_, 
    individuals_);
}

int Search::submitDemes(DB *db) {
  QPair<QString, DBElement*> refMember;
  QHash<QString, QVariant> values;
  QHash<QString, DBElement*> members;

  return ed_.submit(db, refMember, members, values, demes_);
}

int Search::submitParetoFront(DB *db) {
  QPair<QString, DBElement*> refMember;
  QHash<QString, QVariant> values;
  QHash<QString, DBElement*> members;

  return ed_.submit(db, refMember, members, values, paretoFront_);
}

int Search::submitEvolution(DB *db) {
  QPair<QString, DBElement*> refMember;

  QHash<QString, QVariant> values;
  values.insert("Name", name_);
  values.insert("RandSeed", randSeed_);
  values.insert("StartDatetime", startDatetime_.toString(DB::DatetimeFormat));
  values.insert("EndDatetime", endDatetime_.toString(DB::DatetimeFormat));

  QHash<QString, DBElement*> members;

  // Submitting pareto front and old pareto front individuals
  return ed_.submit(db, refMember, members, values, demes_, 
                    paretoFront_ + oldParetoFrontInds_);
}

int Search::submitShallow(DB *db) {
  QHash<QString, QVariant> values;
  values.insert("Name", name_);
  values.insert("RandSeed", randSeed_);
  values.insert("StartDatetime", startDatetime_.toString(DB::DatetimeFormat));
  values.insert("EndDatetime", endDatetime_.toString(DB::DatetimeFormat));

  return ed_.submit(db, values);
}

bool Search::erase() {
  QList<DBElement*> members;

  return ed_.erase(members, searchExperiments_, individuals_, demes_);
}

}