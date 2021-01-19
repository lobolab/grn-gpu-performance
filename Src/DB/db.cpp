// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#include "db.h"
////#include "log.h"

#include <QSqlQuery>
#include <QSqlRecord>
#include <QSqlRelationalTableModel>
#include <QSqlError>
#include <QFile>
#include <QStringList>
#include <QDebug>
#include <QUuid>

namespace LoboLab {

const char *DB::DatetimeFormat = "yyyy-MM-ddThh:mm:ss.zzzZ";

DB::DB()
  : nNestedTrans_(0) {
}

DB::~DB() {
  disconnect();
}

int DB::connect(const QString &fileName, bool inFastMode) {
  int error = 0;

  if (isConnected())
    disconnect();

  if (QFile(fileName).exists()) {
    QString connectionName = QUuid::createUuid().toString();
    db_ = QSqlDatabase::addDatabase("QSQLITE", connectionName);
    db_.setDatabaseName(fileName);
    if (db_.open()) {
      QSqlQuery query(db_);

      if(inFastMode) {
        query.exec("PRAGMA foreign_keys = false;");
        Q_ASSERT_X(query.isActive(), "DB::connect: setting PRAGMA :",
          query.lastError().text().toLatin1());

        query.exec("PRAGMA temp_store = 2;");
        Q_ASSERT_X(query.isActive(), "DB::connect: setting PRAGMA :",
          query.lastError().text().toLatin1());

        query.exec("PRAGMA synchronous = 0;");
        Q_ASSERT_X(query.isActive(), "DB::connect: setting PRAGMA :",
          query.lastError().text().toLatin1());
      } else { // Normal mode
        query.exec("PRAGMA foreign_keys = true;");
        Q_ASSERT_X(query.isActive(), "DB::connect: setting PRAGMA :",
          query.lastError().text().toLatin1());
      }

    } else
      error = 1;
  } else
    error = 2;

  return error;
}

void DB::disconnect() {
  QString connectionName = db_.connectionName();
  db_ = QSqlDatabase();
  QSqlDatabase::removeDatabase(connectionName);
}

bool DB::isConnected() {
  return db_.isOpen();
}

bool DB::createEmptyDB(const QString &fileName) {
  if (db_.isOpen())
    disconnect();

  QFile file(fileName);
  bool ok = file.open(QIODevice::WriteOnly | QIODevice::Truncate);

  Q_ASSERT_X(ok, "DB::createEmptyDB:",
             QString("file error code =%1").arg(file.error()).toLatin1());

  file.close();

  ok &= connect(fileName) == 0;

  return ok;
}

bool DB::beginTransaction() {
  bool ok;

  if (nNestedTrans_ == 0) {
    QSqlQuery query(db_);
    ok = query.exec("BEGIN TRANSACTION;");

    Q_ASSERT_X(ok, "DB::beginTransaction:",
               query.lastError().text().toLatin1());
  } else
    ok = true;

  ++nNestedTrans_;

  return ok;
}

bool DB::endTransaction() {
  bool ok = true;
  --nNestedTrans_;

  if (nNestedTrans_ == 0) {
    QSqlQuery query(db_);
    ok = query.exec("COMMIT;");

    Q_ASSERT_X(ok, "DB::endTransaction:",
               query.lastError().text().toLatin1());
  }

  if (nNestedTrans_ < 0) {
    nNestedTrans_ = 0;
    ok = false;
  }

  return ok;
}

bool DB::rollbackTransaction() {
  bool ok;
  if (nNestedTrans_ > 0) {
    QSqlQuery query(db_);
    ok = query.exec("ROLLBACK;");

    Q_ASSERT_X(ok, "DB::rollbackTransaction:",
               query.lastError().text().toLatin1());

    nNestedTrans_ = 0;
  } else if (nNestedTrans_ < 0) {
    nNestedTrans_ = 0;
    ok = false;
  } else
    ok = true;

  return ok;
}

bool DB::vacuum() {
  bool ok;
  QSqlQuery query("VACUUM", db_);

  ok = query.exec();

  Q_ASSERT_X(ok, QString("DB::vacuum").toLatin1(),
             query.lastError().text().toLatin1());

  return ok;
}

QSqlQuery *DB::newQuery() const {
  QSqlQuery *query = new QSqlQuery(db_);
  query->setForwardOnly(true);

  return query;
}

QSqlQuery *DB::newQuery(const QString &sqlStr) const {
  QSqlQuery *query = new QSqlQuery(sqlStr, db_);
  query->setForwardOnly(true);
  query->exec();

  Q_ASSERT_X(query->isActive(), ("DB::newQuery: " + sqlStr).toLatin1(),
             query->lastError().text().toLatin1());

  return query;
}

QSqlQuery *DB::newQuery(const QString &sqlStr, const QVariant &value) const {
  QSqlQuery *query = new QSqlQuery(NULL, db_);
  query->setForwardOnly(true);
  bool ok = query->prepare(sqlStr);
  query->addBindValue(value);
  ok &= query->exec();

  Q_ASSERT_X(ok, QString("DB::newQuery: %1 value=%2")
             .arg(sqlStr).arg(value.toString()).toLatin1(),
             query->lastError().text().toLatin1());

  return query;
}

QSqlQuery *DB::newQuery(const QString &sqlStr,
                        const QVariantList &values) const {
  QSqlQuery *query = new QSqlQuery(NULL, db_);
  query->setForwardOnly(true);
  bool ok = query->prepare(sqlStr);

  int nValues = values.size();
  for (int i = 0; i<nValues; ++i)
    query->addBindValue(values.at(i));

  ok &= query->exec();

  Q_ASSERT_X(ok, QString("DB::newQuery: %1 nValues=%2")
             .arg(sqlStr).arg(nValues).toLatin1(),
             query->lastError().text().toLatin1());

  return query;
}

QSqlQuery *DB::newTableQuery(const QString &table) const {
  QString sqlStr = "SELECT * FROM " + table;
  QSqlQuery *query = new QSqlQuery(sqlStr, db_);
  query->setForwardOnly(true);
  query->exec();

  Q_ASSERT_X(query->isActive(), ("DB::newQuery: " + sqlStr).toLatin1(),
             query->lastError().text().toLatin1());

  return query;
}

QSqlQuery *DB::newTableQuery(const QString &table, int id) const {
  QString sqlStr = QString("SELECT * FROM " + table + " WHERE Id = %1")
                   .arg(id);
  QSqlQuery *query = new QSqlQuery(sqlStr, db_);
  query->setForwardOnly(true);
  query->exec();

  Q_ASSERT_X(query->isActive(), ("DB::newQuery: " + sqlStr).toLatin1(),
             query->lastError().text().toLatin1());

  return query;
}

QSqlQuery *DB::newTableQuery(const QString &table, const QString &refField,
                             int id) const {
  QString sqlStr = QString("SELECT * FROM " + table + " WHERE " + table + "."
                           + refField + " = %1").arg(id);

  QSqlQuery *query = new QSqlQuery(sqlStr, db_);
  query->setForwardOnly(true);
  query->exec();

  Q_ASSERT_X(query->isActive(), ("DB::newQuery: " + sqlStr).toLatin1(),
             query->lastError().text().toLatin1());

  return query;
}

QSqlQuery *DB::newTableQuery(const QString &table, const QString &refField,
                             int id, const QString &sortField) const {
  QString sqlStr = QString("SELECT * FROM " + table + " WHERE " + table + "."
                           + refField + " = %1 ORDER BY " + sortField).arg(id);

  QSqlQuery *query = new QSqlQuery(sqlStr, db_);
  query->setForwardOnly(true);
  query->exec();

  Q_ASSERT_X(query->isActive(), ("DB::newQuery: " + sqlStr).toLatin1(),
             query->lastError().text().toLatin1());

  return query;
}

QSqlQueryModel *DB::newModel(const QString &sqlStr) const {
  QSqlQuery query(sqlStr, db_);
  query.exec();

  Q_ASSERT_X(query.isActive(), ("DB::newModel: " + sqlStr).toLatin1(),
             query.lastError().text().toLatin1());

  QSqlQueryModel *model = new QSqlQueryModel();
  model->setQuery(query);

  Q_ASSERT_X(!model->lastError().isValid(), QString("DB::newModel: " + sqlStr)
             .toLatin1(), model->lastError().text().toLatin1());

  fetchAllData(model);

  return model;
}

QSqlQueryModel *DB::newModel(const QString &sqlStr, const QVariant &value) const 
{
  QSqlQuery query(NULL, db_);
  bool ok = query.prepare(sqlStr);
  query.addBindValue(value);
  ok &= query.exec();

  Q_ASSERT_X(ok, QString("DB::query: %1 value=%2")
             .arg(sqlStr).arg(value.toString()).toLatin1(),
             query.lastError().text().toLatin1());

  QSqlQueryModel *model = new QSqlQueryModel();
  model->setQuery(query);

  Q_ASSERT_X(!model->lastError().isValid(), QString("DB::newModel: %1")
             .arg(sqlStr).toLatin1(), model->lastError().text().toLatin1());

  fetchAllData(model);

  return model;
}

QSqlQueryModel *DB::newModel(const QString &sqlStr,
                             const QVariantList &values) const {
  QSqlQuery query(NULL, db_);
  bool ok = query.prepare(sqlStr);

  int nValues = values.size();
  for (int i = 0; i<nValues; ++i)
    query.addBindValue(values.at(i));

  ok &= query.exec();

  Q_ASSERT_X(ok, QString("DB::query: %1 nValues=%2")
             .arg(sqlStr).arg(nValues).toLatin1(),
             query.lastError().text().toLatin1());

  QSqlQueryModel *model = new QSqlQueryModel();
  model->setQuery(query);

  Q_ASSERT_X(!model->lastError().isValid(), QString("DB::query: %1")
             .arg(sqlStr).toLatin1(), model->lastError().text().toLatin1());

  fetchAllData(model);

  return model;
}

QSqlQueryModel *DB::newTableModel(const QString &table) const {
  QSqlQueryModel *model = new QSqlQueryModel();
  QString sqlStr = "SELECT * FROM " + table;
  model->setQuery(sqlStr, db_);

  Q_ASSERT_X(!model->lastError().isValid(), QString("DB::newTableModel: " +
             table).toLatin1(), model->lastError().text().toLatin1());

  fetchAllData(model);

  return model;
}

QSqlQueryModel *DB::newTableModel(const QString &table, int id) const {
  QSqlQueryModel *model = new QSqlQueryModel();
  QString sqlStr = QString("SELECT * FROM %1 WHERE %1.Id = %2").arg(table)
                   .arg(id);
  model->setQuery(sqlStr, db_);

  Q_ASSERT_X(!model->lastError().isValid(), QString("DB::newTableModel (id): "
             + table).toLatin1(), model->lastError().text().toLatin1());

  fetchAllData(model);

  return model;
}

QSqlQueryModel *DB::newTableModel(const QString &table, const QString &refField,
                                  int id) const {
  QSqlQueryModel *model = new QSqlQueryModel();
  QString sqlStr = QString("SELECT * FROM %1 WHERE %1.%2 = %3").arg(table).
                   arg(refField).arg(id);
  model->setQuery(sqlStr, db_);

  Q_ASSERT_X(!model->lastError().isValid(), QString("DB::newTableModel (ref):"
             + table).toLatin1(), model->lastError().text().toLatin1());

  fetchAllData(model);

  return model;
}

QSqlQueryModel *DB::newTableModel(const QString &table, const QString &refField,
                                  int id, const QString &sortField) const {
  QSqlQueryModel *model = new QSqlQueryModel();
  QString sqlStr = QString("SELECT * FROM %1 WHERE %1.%2 = %3 ORDER BY %4")
                   .arg(table).arg(refField).arg(id).arg(sortField);
  model->setQuery(sqlStr, db_);

  Q_ASSERT_X(!model->lastError().isValid(), QString("DB::newTableModel (ref):"
             + table).toLatin1(), model->lastError().text().toLatin1());

  fetchAllData(model);

  return model;
}

QSqlQueryModel *DB::newTableModel(const QString &table,
                                  const QString &sortField) const {
  QSqlQueryModel *model = new QSqlQueryModel();
  QString sqlStr = "SELECT * FROM " + table + " ORDER BY " + sortField;
  model->setQuery(sqlStr, db_);

  Q_ASSERT_X(!model->lastError().isValid(), QString("DB::newTableModel: " +
             table).toLatin1(), model->lastError().text().toLatin1());

  fetchAllData(model);

  return model;
}

QSqlQueryModel *DB::newTableModel(const QString &table,
                                  const QMultiHash<QString, QStringList> &relations) const {
  QSqlQueryModel *model = new QSqlQueryModel();

  QString sqlStr = "SELECT " + table + ".Id," + table + ".Name";

  for (QMultiHash<QString, QStringList>::const_iterator i = relations.constBegin();
       i != relations.constEnd(); ++i) {
    sqlStr += ", " + i.value().at(2) + "." + i.value().at(1) + " AS '" + i.value().at(0) + "'";
  }

  sqlStr += " FROM " + table;

  for (QMultiHash<QString, QStringList>::const_iterator i = relations.constBegin();
       i != relations.constEnd(); ++i) {
    sqlStr += " INNER JOIN " + i.key() + " " + i.value().at(2) + " ON " + table + "." +
      i.value().at(0) + " = " + i.value().at(2) + ".Id";
  } 
  
  model->setQuery(sqlStr, db_);

  Q_ASSERT_X(!model->lastError().isValid(), QString("DB::newTableModel "
             "(relations): " + table).toLatin1(), model->lastError().text()
             .toLatin1());

  fetchAllData(model);

  return model;
}

void DB::updateModel(QSqlQueryModel *model) const { 
  model->setQuery(model->query().lastQuery(), db_);
  fetchAllData(model);
}

int DB::getModelIdRow(QSqlQueryModel *model, int id) const {
  int row;
  QModelIndex start = model->index(0, 0);
  QModelIndexList indexList = model->match(start, Qt::DisplayRole, id, 1,
                              Qt::MatchExactly | Qt::MatchWrap);

  Q_ASSERT(indexList.count() <= 1);

  if (indexList.count() == 0)
    row = -1;
  else
    row = indexList.first().row();

  return row;
}

int DB::insertRow(const QString &table, const QHash<QString, QVariant> &values,
                  bool ignore) {
  bool ok;
  QSqlQuery query(NULL, db_);

  if (!values.isEmpty()) {
    QString columnsStr(" (");
    QString valuesStr(" (");
    for (QHash<QString, QVariant>::const_iterator i = values.constBegin();
         i!=values.constEnd(); ++i) {
      columnsStr.append(i.key() + ", ");
      valuesStr.append("?, ");
    }
    columnsStr.replace(columnsStr.length()-2, 1, ')');
    valuesStr.replace(valuesStr.length()-2, 1, ')');

    QString insertStr;
    if (ignore)
      insertStr = "INSERT OR IGNORE INTO ";
    else
      insertStr = "INSERT INTO ";

    ok = query.prepare(insertStr + table + columnsStr + "VALUES" +
                       valuesStr);

    Q_ASSERT_X(ok, QString("DB::insertRow: %1")
               .arg(table).toLatin1(), query.lastError().text().toLatin1());

    for (QHash<QString, QVariant>::const_iterator i = values.constBegin();
         i!=values.constEnd(); ++i)
      query.addBindValue(i.value());

    ok &= query.exec();
  } else
    ok = query.exec("INSERT INTO " + table + " DEFAULT VALUES");

  Q_ASSERT_X(ok, QString("DB::insertRow: %1")
             .arg(table).toLatin1(), query.lastError().text().toLatin1());

  return query.lastInsertId().toInt();
}

bool DB::updateRow(const QString &table, int id,
                   const QHash<QString, QVariant> &values) const {
  bool ok;
  QSqlQuery query(NULL, db_);

  if (!values.isEmpty()) {
    QString sql = "UPDATE " + table + " SET ";
    for (QHash<QString, QVariant>::const_iterator i = values.constBegin();
         i!=values.constEnd(); ++i)
      sql += i.key() + "=?, ";

    sql.remove(sql.length()-2, 2);
    sql += QString(" WHERE Id=%1").arg(id);

    ok = query.prepare(sql);

    for (QHash<QString, QVariant>::const_iterator i = values.constBegin();
         i!=values.constEnd(); ++i)
      query.addBindValue(i.value());

    ok &= query.exec();
  } else
    ok = true;

  Q_ASSERT_X(ok, QString("DB::updateRow: %1")
             .arg(table).toLatin1(), query.lastError().text().toLatin1());

  return ok;
}

bool DB::removeRow(const QString &table, int id) {
  bool ok;
  QString sql = QString("DELETE FROM %1 WHERE Id=%2").arg(table).arg(id);
  QSqlQuery query(sql, db_);

  ok = query.exec();

  Q_ASSERT_X(ok, QString("DB::removeRow: %1")
             .arg(table).toLatin1(), query.lastError().text().toLatin1());

  return ok;
}

void DB::fetchAllData(QSqlQueryModel *model) const {
  QModelIndex invalidIndex;
  while (model->canFetchMore(invalidIndex))
    model->fetchMore(invalidIndex);
}

bool DB::exist(const QString &table) const {
  QSqlQuery query(NULL, db_);
  bool ok = query.exec(QString("SELECT name FROM sqlite_master WHERE "
                               "type='table' AND name='%1'").arg(table));
  ok &= query.exec();

  Q_ASSERT_X(ok, QString("DB::exist table=%1").arg(table).toLatin1(),
             query.lastError().text().toLatin1());

  return query.first();
}

bool DB::exist(const QString &table, const QString &field,
               const QVariant &content) const {
  QSqlQuery query(NULL, db_);
  bool ok;

  if (content.isNull())
    ok = query.exec(QString("SELECT COUNT(*) FROM %1 WHERE %2 IS null")
                    .arg(table).arg(field));
  else {
    ok = query.prepare(QString("SELECT COUNT(%2) FROM %1 WHERE %2=?")
                       .arg(table).arg(field));
    query.addBindValue(content);
    ok &= query.exec();
  }

  Q_ASSERT_X(ok, QString("DB::exist table=%1 field=%2 content=%3")
             .arg(table).arg(field).arg(content.toString()).toLatin1(),
             query.lastError().text().toLatin1());

  ok = query.first();
  Q_ASSERT_X(ok, QString("DB::exist query.first() table=%1 field=%2 content=%3")
             .arg(table).arg(field).arg(content.toString()).toLatin1(),
             query.lastError().text().toLatin1());

  int nRows = query.value(0).toInt();

  return nRows > 0;
}

int DB::existId(const QString &table, const QString &field,
                const QVariant &content) const {
  QSqlQuery query(NULL, db_);

  bool ok;

  if (content.isNull())
    ok = query.exec(QString("SELECT %1.Id FROM %1 WHERE %2 IS null")
                    .arg(table).arg(field));
  else {
    ok = query.prepare(QString("SELECT %1.Id FROM %1 WHERE %2=?")
                       .arg(table).arg(field));
    query.addBindValue(content);
    ok &= query.exec();
  }

  Q_ASSERT_X(ok, QString("DB::existId table=%1 field=%2 content=%3")
             .arg(table).arg(field).arg(content.toString()).toLatin1(),
             query.lastError().text().toLatin1());

  int id;
  if (query.next())
    id = query.value(0).toInt();
  else
    id = 0;

  Q_ASSERT(!query.next());

  return id;
}

int DB::getNumRows(const QString &table) const {
  QSqlQuery query(QString("SELECT COUNT(*) FROM %1").arg(table), db_);
  query.exec();

  Q_ASSERT_X(query.isActive(), ("DB::nRows: " + table).toLatin1(),
             query.lastError().text().toLatin1());

  query.next();
  return query.value(0).toInt();
}

//void DB::checkError(QAbstractItemModel *model)
//{
//	qDebug() << ((QSqlTableModel*)model)->lastError();
//}

void DB::emptyTables(const QList<QString> &tables) const {
  bool ok = true;
  QSqlQuery query(db_);
  for (int i = 0; i<tables.size(); ++i) {
    ok &= query.exec("DELETE FROM " + tables.at(i));

    Q_ASSERT_X(ok, QString("DB::emptyTables: table=%1")
               .arg(tables.at(i)).toLatin1(), query.lastError().text().toLatin1());
  }

  Q_ASSERT_X(ok, "DB::emptyTables", query.lastError().text().toLatin1());
}

void DB::emptyTablesNoKeys(const QList<QString> &tables) const {
  QSqlQuery query(db_);
  bool ok = query.exec("PRAGMA foreign_keys = false;");
  Q_ASSERT_X(ok, "DB::emptyTablesNoKeys: removing foreing_keys pragma",
             query.lastError().text().toLatin1());

  for (int i = 0; i<tables.size(); ++i) {
    ok &= query.exec("DELETE FROM " + tables.at(i));

    Q_ASSERT_X(ok, QString("DB::emptyTablesNoKeys: table=%1")
               .arg(tables.at(i)).toLatin1(), query.lastError().text().toLatin1());
  }

  ok &= query.exec("PRAGMA foreign_keys = true;");
  Q_ASSERT_X(ok, "DB::emptyTablesNoKeys: restoring foreing_keys pragma",
             query.lastError().text().toLatin1());
}

}