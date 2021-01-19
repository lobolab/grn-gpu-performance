// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#pragma once

namespace LoboLab {

class DB;

// Persistence class for hierarchical elements stored in the DB
class DBElement {
  friend class DBElementData;

 protected:
  virtual ~DBElement() {}

  virtual int id() const = 0;
  virtual int submit(DB *db) = 0;
  virtual bool erase() = 0;
};

} // namespace LoboLab
