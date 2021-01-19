// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#pragma once

#include "DB/dbelementdata.h"

namespace LoboLab {

class Search;
class Experiment;

class SearchExperiment : public DBElement {
  friend class Search;

 public:
  inline Search *search() const { return search_; };
  inline Experiment *experiment() const { return experiment_; }
  inline void setExperiment(Experiment *experiment) { 
    experiment_ = experiment; 
  }

 protected:
  inline virtual int id() const { return ed_.id(); };
  virtual int submit(DB *db);
  virtual bool erase();

 private:
  SearchExperiment(Search *s, Experiment *e);
  SearchExperiment(Search *s, const DBElementData &ref);
  ~SearchExperiment();

  SearchExperiment(const SearchExperiment &source, Search *s = NULL,
                   bool maintainId = true);
  SearchExperiment &operator=(const SearchExperiment &source);

  void load();

  Search *search_;
  Experiment *experiment_;

  DBElementData ed_;

// Persistence fields
 public:
  enum {
    FSearch = 1,
    FExperiment
  };
};

} // namespace LoboLab
