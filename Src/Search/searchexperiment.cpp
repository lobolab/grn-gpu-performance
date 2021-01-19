// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#include "searchexperiment.h"
#include "search.h"
#include "Experiment/experiment.h"

namespace LoboLab {

SearchExperiment::SearchExperiment(Search *s, Experiment *e)
  : search_(s), experiment_(e), ed_("SearchExperiment") {
}

SearchExperiment::SearchExperiment(Search *s, const DBElementData &ref)
  : search_(s), ed_("SearchExperiment", ref) {
  load();
}

SearchExperiment::~SearchExperiment() {
  delete experiment_;
}

SearchExperiment::SearchExperiment(const SearchExperiment &source, Search *s,
                                   bool maintainId)
  : search_(s), ed_(source.ed_, maintainId) {
  experiment_ = new Experiment(*source.experiment_);
}


// Persistence methods

void SearchExperiment::load() {
  experiment_ = new Experiment(ed_.loadValue(FExperiment).toInt(), ed_.db());

  ed_.loadFinished();
}

int SearchExperiment::submit(DB *db) {
  QPair<QString, DBElement*> refMember("Search", search_);

  QHash<QString, DBElement*> members;
  members.insert("Experiment", experiment_);

  return ed_.submit(db, refMember, members);
}

bool SearchExperiment::erase() {
  return ed_.erase();
}

}