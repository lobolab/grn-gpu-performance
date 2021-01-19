// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#pragma once

#include "DB/db.h"
#include <QObject>

namespace LoboLab {
class DB;
class Search;

class MainCmd : public QObject {
  Q_OBJECT

 public:
  MainCmd(QObject *parent = NULL);
  ~MainCmd(void);

 private slots:
  void run();

 private:
  void quit(int ret);
  int runSearch();
  void closeDB();
  bool connectDB(DB &db, const QString &dbFileName);

  DB db_;
  Search *search_;
};

} // namespace LoboLab
