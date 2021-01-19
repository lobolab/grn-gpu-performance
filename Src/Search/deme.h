// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#pragma once

#include "DB/dbelementdata.h"
#include <QDateTime>

namespace LoboLab {

class Search;
class Generation;
class Individual;

class Deme : public DBElement {
  friend class Search;
 public:
  ~Deme();

  inline Search *search() const {return search_;}  
  inline const QList<Generation*> &generations() const { return generations_; }
  inline QList<Generation*> &generations() { return generations_; }
  
  Generation *createNextGeneration();
 
  inline virtual int id() const { return ed_.id(); };
  virtual int submit(DB *db);
  int submitShallow(DB *db);
  virtual bool erase();

 private:
  explicit Deme(Search *s);
  Deme(Search *s, const QHash<int, Individual*> &individualsIdMap, 
       const DBElementData &ref);
  Deme(const Deme &source, Search *s = NULL, bool maintainId = true);
  Deme &operator=(const Deme &source);

  void load(const QHash<int, Individual*> &individualsIdMap);
  void loadGenerations(const QHash<int, Individual*> &individualsIdMap);

  Search *search_;
  QList<Generation*> generations_;

  DBElementData ed_;

// Persistence fields
 public:
  enum {
    FSearch = 1
  };
};

} // namespace LoboLab
