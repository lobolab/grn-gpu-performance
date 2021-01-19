// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#pragma once

#include <QSqlQuery>
#include <QHash>
#include <QVariant>
#include <QPair>

#include "dbelement.h"

namespace LoboLab {

class DB;
class DBExp;
class DBMod;
class DBSea;

// Persistence class for elements stored in the DB
class DBElementData {
 public:
  DBElementData(const char *elementName);
  DBElementData(const char *elementName, int id, DB *db);
  DBElementData(const char *elementName, const DBElementData &ref);
  virtual ~DBElementData();

  DBElementData(const DBElementData &source, bool maintainId = true);
  virtual DBElementData &operator=(const DBElementData &source);

  QVariant loadValue(int index);
  void loadReferences(const QString &refName);
  void loadReferences(const QString &refName, const QString &sortField);
  void loadReferences(const QString &refName, const QString &refField,
                      const QString &referencedName);
  void loadReferencesIndirect(const QString &refName1, const QString &refName2, 
                              const QString &refName3);
  void loadReferencesIndirect(const QString &refName1, const QString &refName2, 
                              const QString &refName3, const QString &refName4);
  bool nextReference();
  QVariant loadRefValue(int index);
  void loadExtElement(const QString &elementName, int id);
  QVariant loadExtValue(int index);
  void loadFinished();

  int submit(DB *db, const QHash<QString, QVariant> &values);

  int submit(DB *db, const QHash<QString, DBElement*> &members,
             const QHash<QString, QVariant> &values = (QHash<QString, QVariant>()));

  int submit(DB *db, const QPair<QString, DBElement*> &refMember,
             const QHash<QString, QVariant> &values = (QHash<QString, QVariant>()));

  int submit(DB *db, const QPair<QString, DBElement*> &refMember,
             const QHash<QString, DBElement*> &members,
             const QHash<QString, QVariant> &values = (QHash<QString, QVariant>()));

  
  template <class T>
  int submit(DB *dbin, const QHash<QString, QVariant> &values,
             const QList<T*> &references);

  template <class T>
  int submit(DB *db, const QPair<QString, DBElement*> &refMember,
             const QHash<QString, DBElement*> &members,
             const QHash<QString, QVariant> &values,
             const QList<T*> &references);

  template <class T, class U>
  int submit(DB *db, const QPair<QString, DBElement*> &refMember,
             const QHash<QString, DBElement*> &members,
             const QHash<QString, QVariant> &values,
             const QList<T*> &references1,
             const QList<U*> &references2);

  template <class T, class U, class V>
  int submit(DB *db, const QPair<QString, DBElement*> &refMember,
             const QHash<QString, DBElement*> &members,
             const QHash<QString, QVariant> &values,
             const QList<T*> &references1,
             const QList<U*> &references2,
             const QList<V*> &references3);

  virtual bool erase();
  bool erase(const QList<DBElement*> &members);

  template <class T>
  int erase(const QList<DBElement*> &members,
            const QList<T*> &references);

  template <class T, class U>
  int erase(const QList<DBElement*> &members,
            const QList<T*> &references1,
            const QList<U*> &references2);

  template <class T, class U, class V>
  int erase(const QList<DBElement*> &members,
            const QList<T*> &references1,
            const QList<U*> &references2,
            const QList<V*> &references3);

  void removeReference(DBElement *reference);

  inline int id() const {
    return id_;
  }
  inline DB *db() const {
    return db_;
  }

 private:
  bool startSubmit(DB *dbin);
  void startErase();
  bool submit(const QHash<QString, QVariant> &values);
  bool submit(const QHash<QString, DBElement*> &members,
              const QHash<QString, QVariant> &values);
  bool submit(const QPair<QString, DBElement*> &refMember,
              const QHash<QString, QVariant> &values);
  bool submit(const QPair<QString, DBElement*> &refMember,
              const QHash<QString, DBElement*> &members,
              const QHash<QString, QVariant> &values);
  template <class T>
  bool submit(const QList<T*> &references);

  void endSubmit(bool ok);
  void endErase(bool ok);

  bool eraseElement();
  bool eraseMembers(const QList<DBElement*> &members);
  template <class T>
  bool eraseReferences(const QList<T*> &references);
  bool eraseRemovedReferences();
  void deleteRemovedReferences();

  const char *elementName_; // Using char* avoids the overhead of creating a
                           // QString, especially for those classes that will 
                           // not use the db.
  int id_;

  DB *db_;
  QSqlQuery *query_;
  QSqlQuery *queryReferences_;
  QSqlQuery *queryExt_;

  QList<DBElement*> removedReferences_;
};


// The template functions definitions are here because of the following:
// http://www.parashift.com/c++-faq-lite/templates.html#faq-35.12

template <class T>
int DBElementData::submit(DB *dbin, const QHash<QString, QVariant> &values,
                          const QList<T*> &references) {
  bool ok = startSubmit(dbin);
  ok &= submit(values);
  ok &= submit(references);
  endSubmit(ok);

  return id_;
}

template <class T>
int DBElementData::submit(DB *dbin, const QPair<QString, DBElement*> &refMember,
                          const QHash<QString, DBElement*> &members,
                          const QHash<QString, QVariant> &values,
                          const QList<T*> &references) {
  bool ok = startSubmit(dbin);
  ok &= submit(refMember, members, values);
  ok &= submit(references);
  endSubmit(ok);

  return id_;
}

template <class T, class U>
int DBElementData::submit(DB *dbin, const QPair<QString, DBElement*> &refMember,
                          const QHash<QString, DBElement*> &members,
                          const QHash<QString, QVariant> &values,
                          const QList<T*> &references1,
                          const QList<U*> &references2) {
  bool ok = startSubmit(dbin);
  ok &= submit(refMember, members, values);
  ok &= submit(references1);
  ok &= submit(references2);
  endSubmit(ok);

  return id_;
}

template <class T, class U, class V>
int DBElementData::submit(DB *dbin, const QPair<QString, DBElement*> &refMember,
                          const QHash<QString, DBElement*> &members,
                          const QHash<QString, QVariant> &values,
                          const QList<T*> &references1,
                          const QList<U*> &references2,
                          const QList<V*> &references3) {
  bool ok = startSubmit(dbin);
  ok &= submit(refMember, members, values);
  ok &= submit(references1);
  ok &= submit(references2);
  ok &= submit(references3);
  endSubmit(ok);

  return id_;
}

template <class T>
bool DBElementData::submit(const QList<T*> &references) {
  bool ok = true;

  int n = references.size();
  for (int i = 0; i < n; ++i)
    ok &= ((DBElement*)(references.at(i)))->submit(db_) > 0;

  return ok;
}

template <class T>
int DBElementData::erase(const QList<DBElement*> &members,
                         const QList<T*> &references) {
  startErase();
  bool ok = eraseRemovedReferences();
  ok &= eraseReferences(references);
  ok &= eraseElement();
  ok &= eraseMembers(members);
  endErase(ok);

  return ok;
}

template <class T, class U>
int DBElementData::erase(const QList<DBElement*> &members,
                         const QList<T*> &references1,
                         const QList<U*> &references2) {
  startErase();
  bool ok = eraseRemovedReferences();
  ok &= eraseReferences(references1);
  ok &= eraseReferences(references2);
  ok &= eraseElement();
  ok &= eraseMembers(members);
  endErase(ok);

  return ok;
}

template <class T, class U, class V>
int DBElementData::erase(const QList<DBElement*> &members,
                         const QList<T*> &references1,
                         const QList<U*> &references2,
                         const QList<V*> &references3) {
  startErase();
  bool ok = eraseRemovedReferences();
  ok &= eraseReferences(references1);
  ok &= eraseReferences(references2);
  ok &= eraseReferences(references3);
  ok &= eraseElement();
  ok &= eraseMembers(members);
  endErase(ok);

  return ok;
}

template <class T>
bool DBElementData::eraseReferences(const QList<T*> &references) {
  bool ok = true;
  int n = references.size();
  for (int i = 0; i < n; ++i)
    ok &= ((DBElement*)(references.at(i)))->erase();

  return ok;
}

} // namespace LoboLab
