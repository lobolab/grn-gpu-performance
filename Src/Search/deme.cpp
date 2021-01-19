// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#include "deme.h"
#include "search.h"
#include "searchparams.h"
#include "generation.h"
#include "individual.h"
#include "DB/db.h"

namespace LoboLab {

Deme::Deme(Search *s)
  : search_(s), ed_("Deme") {
    generations_.reserve(search_->searchParams()->nGenerations);
}

Deme::Deme(Search *s, const QHash<int, Individual*> &individualsIdMap, 
           const DBElementData &ref)
    : search_(s), ed_("Deme", ref) {
  load(individualsIdMap);
}

Deme::~Deme() {
  int n = generations_.size();
  for (int i = 0; i < n; ++i)
    delete generations_.at(i);
}

Generation *Deme::createNextGeneration() {
  int ind;
  if (generations_.isEmpty())
    ind = 0;
  else
    ind = generations_.last()->ind() + 1;

  Generation *gen = new Generation(this, ind);
  generations_.append(gen);
  return gen;
}

// Persistence methods

void Deme::load(const QHash<int, Individual*> &individualsIdMap) {
  loadGenerations(individualsIdMap);

  ed_.loadFinished();
}

void Deme::loadGenerations(const QHash<int, Individual*> &individualsIdMap) {
  ed_.loadReferences("Generation");
  while (ed_.nextReference()) {
    Generation *ele = new Generation(this, individualsIdMap, ed_);
    generations_.append(ele);
  }
}

int Deme::submit(DB *db) {
  QPair<QString, DBElement*> refMember("Search", search_);

  QHash<QString, DBElement*> members;
  QHash<QString, QVariant> values;

  return ed_.submit(db, refMember, members, values, generations_);
}

int Deme::submitShallow(DB *db) {
  QPair<QString, DBElement*> refMember("Search", search_);
  
  return ed_.submit(db, refMember);
}

bool Deme::erase() {
  QList<DBElement*> members;

  return ed_.erase(members, generations_);
}

}