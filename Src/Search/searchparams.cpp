// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#include "searchparams.h"
#include "DB/db.h"

namespace LoboLab {

SearchParams::SearchParams(int id, DB *db)
  : ed("SearchParams", id, db) {
  load();
}

SearchParams::~SearchParams() {
}

SearchParams::SearchParams(const SearchParams &source, bool maintainId)
  : ed(source.ed, maintainId) {
  copy(source);
}

SearchParams &SearchParams::operator=(const SearchParams &source) {
  ed = source.ed;
  copy(source);

  return *this;
}

void SearchParams::copy(const SearchParams &source) {
  name = source.name;
  nDemes = source.nDemes;
  demesSize = source.demesSize;
  nGenerations = source.nGenerations;
  migrationPeriod = source.migrationPeriod;
  saveIndividuals = source.saveIndividuals;
  testModel = source.testModel;
}

// Persistence methods

void SearchParams::load() {
  name = ed.loadValue(FName).toString();
  nDemes = ed.loadValue(FNumDemes).toInt();
  demesSize = ed.loadValue(FDemesSize).toInt();
  nGenerations = ed.loadValue(FNumGenerations).toInt();
  migrationPeriod = ed.loadValue(FMigrationPeriod).toInt();
  saveIndividuals = ed.loadValue(FSaveIndividuals).toInt();
  testModel = ed.loadValue(FTestModel).toBool();

  ed.loadFinished();
}

int SearchParams::submit(DB *db) {
  QHash<QString, QVariant> values;
  values.insert("Name", name);
  values.insert("NumDemes", nDemes);
  values.insert("DemesSize", demesSize);
  values.insert("MigrationPeriod", migrationPeriod);
  values.insert("SaveIndividuals", saveIndividuals);
  values.insert("TestModel", testModel);

  return ed.submit(db, values);
}

bool SearchParams::erase() {
  return ed.erase();
}
}
