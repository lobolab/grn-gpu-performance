// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#include "dbelementdata.h"
#include "dbelement.h"
#include "db.h"
#include "Common/log.h"

namespace LoboLab {

DBElementData::DBElementData(const char *name)
  : elementName_(name), id_(0), db_(NULL), query_(NULL), queryReferences_(NULL),
    queryExt_(NULL) {
}

DBElementData::DBElementData(const char*name, int eleId, DB *d)
  : elementName_(name), id_(eleId), db_(d), query_(NULL), queryReferences_(NULL),
    queryExt_(NULL) {
}

// This avoids a duplicate access to the database
DBElementData::DBElementData(const char*name, const DBElementData &ref)
  : elementName_(name), id_(ref.queryReferences_->value(0).toInt()),
    db_(ref.db_), query_(new QSqlQuery(*(ref.queryReferences_))),
    queryReferences_(NULL), queryExt_(NULL) {
}

DBElementData::~DBElementData() {
  Q_ASSERT(!query_ && !queryReferences_ && !queryExt_);

  loadFinished();
  deleteRemovedReferences();
}

DBElementData::DBElementData(const DBElementData &source, bool maintainId)
  : elementName_(source.elementName_), db_(source.db_), query_(NULL),
    queryReferences_(NULL), queryExt_(NULL) {
  if (maintainId)
    id_ = source.id_;
  else
    id_ = 0;
}

DBElementData &DBElementData::operator=(const DBElementData &source) {
  loadFinished();
  id_ = source.id_;
  db_ = source.db_;

  return *this;
}


void DBElementData::loadFinished() {
  delete query_;
  query_ = NULL;
  delete queryReferences_;
  queryReferences_ = NULL;
  delete queryExt_;
  queryExt_ = NULL;
}

QVariant DBElementData::loadValue(int index) {
  Q_ASSERT(id_);

  if (id_) {
    if (!query_) {
      query_ = db_->newTableQuery(elementName_, id_);
      query_->next();
    }

    return query_->value(index);
  } else
    return QVariant();
}

void DBElementData::loadReferences(const QString &refName) {
  if (queryReferences_)
    delete queryReferences_;

  queryReferences_ = db_->newTableQuery(refName, elementName_, id_);
}

void DBElementData::loadReferences(const QString &refName,
                                   const QString &sortField) {
  if (queryReferences_)
    delete queryReferences_;

  queryReferences_ = db_->newTableQuery(refName, elementName_, id_, sortField);
}

void DBElementData::loadReferences(const QString &refName,
                                   const QString &refField, const QString &referencedName) {
  if (queryReferences_)
    delete queryReferences_;

  QString sql = QString("SELECT DISTINCT %1.* FROM %1 "
                        "INNER JOIN %3 ON %1.%2 = %3.Id WHERE %3.%4 = %5;").arg(refName)
                .arg(refField).arg(referencedName).arg(elementName_).arg(id_);

  queryReferences_ = db_->newQuery(sql);
}

void DBElementData::loadReferencesIndirect(const QString &refName1,
    const QString &refName2, const QString &refName3) {
  if (queryReferences_)
    delete queryReferences_;

  QString sql = QString("SELECT DISTINCT %1.* FROM %1 "
                        "INNER JOIN %2 ON %2.%1 = %1.Id "
                        "INNER JOIN %3 ON %3.Id = %2.%3 "
                        "WHERE %3.%4 = %5;").arg(refName1).arg(refName2).
                        arg(refName3).arg(elementName_).arg(id_);

  queryReferences_ = db_->newQuery(sql);
}

void DBElementData::loadReferencesIndirect(const QString &refName1,
    const QString &refName2, const QString &refName3, const QString &refName4) {
  if (queryReferences_)
    delete queryReferences_;

  QString sql = QString("SELECT DISTINCT %1.* FROM %1 "
                        "INNER JOIN %2 ON %2.%1 = %1.Id "
                        "INNER JOIN %3 ON %3.Id = %2.%3 "
                        "INNER JOIN %4 ON %4.Id = %3.%4 "
                        "WHERE %4.%5 = %6;").arg(refName1).arg(refName2).
                        arg(refName3).arg(refName4).arg(elementName_).arg(id_);

  queryReferences_ = db_->newQuery(sql);
}

bool DBElementData::nextReference() {
  if (queryReferences_)
    return queryReferences_->next();
  else
    return 0;
}

QVariant DBElementData::loadRefValue(int index) {
  Q_ASSERT(queryReferences_);

  if (queryReferences_)
    return queryReferences_->value(index);
  else
    return QVariant();
}

void DBElementData::loadExtElement(const QString &elementName, int id) {
  if (queryExt_)
    delete queryExt_;

  queryExt_ = db_->newTableQuery(elementName, id);
  queryExt_->next();
}

QVariant DBElementData::loadExtValue(int index) {
  Q_ASSERT(queryExt_);

  if (queryExt_)
    return queryExt_->value(index);
  else
    return QVariant();
}

int DBElementData::submit(DB *dbin, const QHash<QString, QVariant> &values) {
  bool ok = startSubmit(dbin);
  ok &= submit(values);
  endSubmit(ok);

  return id_;
}

int DBElementData::submit(DB *dbin, const QHash<QString, DBElement*> &members,
                          const QHash<QString, QVariant> &values) {
  bool ok = startSubmit(dbin);
  ok &= submit(members, values);
  endSubmit(ok);

  return id_;
}

int DBElementData::submit(DB *dbin, const QPair<QString, DBElement*> &refMember,
                          const QHash<QString, QVariant> &values) {
  bool ok = startSubmit(dbin);
  ok &= submit(refMember, values);
  endSubmit(ok);

  return id_;
}

int DBElementData::submit(DB *dbin, const QPair<QString, DBElement*> &refMember,
                          const QHash<QString, DBElement*> &members,
                          const QHash<QString, QVariant> &values) {
  bool ok = startSubmit(dbin);
  ok &= submit(refMember, members, values);
  endSubmit(ok);

  return id_;
}

bool DBElementData::submit(const QHash<QString, QVariant> &values) {
  bool ok;

  if (id_>0)
    ok = db_->updateRow(elementName_, id_, values);
  else {
    id_ = db_->insertRow(elementName_, values);
    ok = id_ > 0;
  }

  return ok;
}

bool DBElementData::submit(const QHash<QString, DBElement*> &members,
                           const QHash<QString, QVariant> &values) {
  QHash<QString, QVariant> dbValues = values;

  for (QHash<QString, DBElement*>::const_iterator i = members.constBegin();
       i != members.constEnd(); ++i)
    dbValues.insert(i.key(), i.value()->submit(db_));

  return submit(dbValues);
}

bool DBElementData::submit(const QPair<QString, DBElement*> &refMember,
                           const QHash<QString, QVariant> &values) {
  QHash<QString, QVariant> dbValues = values;

  if (!refMember.first.isEmpty())
    dbValues.insert(refMember.first, refMember.second->id());

  return submit(dbValues);
}

bool DBElementData::submit(const QPair<QString, DBElement*> &refMember,
                           const QHash<QString, DBElement*> &members,
                           const QHash<QString, QVariant> &values) {
  QHash<QString, QVariant> dbValues = values;

  if (!refMember.first.isEmpty())
    dbValues.insert(refMember.first, refMember.second->id());

  for (QHash<QString, DBElement*>::const_iterator i = members.constBegin();
       i != members.constEnd(); ++i)
    dbValues.insert(i.key(), i.value()->submit(db_));

  return submit(dbValues);
}

bool DBElementData::startSubmit(DB *dbin) {
  if (db_ != dbin) {
    db_ = dbin;
    id_ = 0;
  }

  Q_ASSERT(db_);

  bool ok;
  if (db_) {
    db_->beginTransaction();
    ok = true;
  } else
    ok = false;

  ok &= eraseRemovedReferences();

  return ok;
}

void DBElementData::startErase() {
  if (db_)
    db_->beginTransaction();
}

void DBElementData::endSubmit(bool ok) {
  if (ok)
    ok &= db_->endTransaction();
  else {
    Log::write() << "DBElementData:endSubmit: error: " <<
                 db_->lastError().text() << endl;
    db_->rollbackTransaction();
    Q_ASSERT(false);
  }
}

void DBElementData::endErase(bool ok) {
  if (db_) {
    if (ok)
      ok &= db_->endTransaction();
    else {
      Log::write() << "DBElementData:endErase: error" <<
                   db_->lastError().text() << endl;
      db_->rollbackTransaction();
      Q_ASSERT(false);
    }
  }
}

bool DBElementData::erase() {
  bool ok = eraseRemovedReferences();
  ok &= eraseElement();
  return ok;
}

bool DBElementData::erase(const QList<DBElement*> &members) {
  bool ok = eraseRemovedReferences();
  ok &= eraseElement();
  ok &= eraseMembers(members);
  return ok;
}

bool DBElementData::eraseElement() {
  bool ok;

  if (id_>0) {
    ok = db_->removeRow(elementName_, id_);
    id_ = 0;
  } else
    ok = true;

  return ok;
}

bool DBElementData::eraseMembers(const QList<DBElement*> &members) {
  bool ok = true;
  int n = members.size();
  for (int i = 0; i < n; ++i)
    if (DBElement *ele = members.at(i))
      ok &= ele->erase();

  return ok;
}

void DBElementData::removeReference(DBElement *reference) {
  removedReferences_.append(reference);
}

bool DBElementData::eraseRemovedReferences() {
  bool ok = true;
  int n = removedReferences_.size();
  for (int i = 0; i < n; ++i) {
    DBElement *member = removedReferences_.at(i);
    ok &= member->erase();
    delete member;
  }

  removedReferences_.clear();

  Q_ASSERT(ok);
  return ok;
}

void DBElementData::deleteRemovedReferences() {
  int n = removedReferences_.size();
  for (int i = 0; i < n; ++i)
    delete removedReferences_.at(i);

  removedReferences_.clear();
}

}