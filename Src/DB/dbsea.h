// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#pragma once

namespace LoboLab {

class DB;

class DBSea {
 public:
  static bool buildDB(DB *db);

 private:
  DBSea();
  virtual ~DBSea();

};

} // namespace LoboLab
