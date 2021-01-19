// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#include "dbsea.h"

#include <QDebug>

namespace LoboLab {

DBSea::DBSea() {
}

DBSea::~DBSea() {
}

bool DBSea::buildDB(DB *db) {
  // TODO
  Q_ASSERT(false);
  return false;


  //char *dump = "";
  //QStringList sqlSts = QString(dump).split(";", QString::SkipEmptyParts);
  //
  //QSqlQuery query(db);
  //int nSts = sqlSts.size();
  //bool ok = true;
  //for(int i=0; i<nSts; ++i)
  //{
  //	ok &= query.exec(sqlSts.at(i));
  //	Q_ASSERT_X(ok, QString("DBMod::buildDBMod: %1").arg(sqlSts.at(i)).toLatin1(),
  //		query.lastError().text().toLatin1());
  //}

  //return ok;
}

}