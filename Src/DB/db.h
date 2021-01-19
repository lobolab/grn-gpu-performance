// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#pragma once

#include <QString>
#include <QHash>
#include <QVariant>
#include <QSqlQuery>
#include <QSqlQueryModel>
#include <QSqlDatabase>
#include <QSqlError>

namespace LoboLab {

class DB {
 public:
  DB();
  virtual ~DB();

  int connect(const QString &fileName, bool inFastMode = false);
  void disconnect();
  bool isConnected();
  bool createEmptyDB(const QString &fileName);
  bool vacuum();

  inline QString fileName() {
    return db_.databaseName();
  };
  inline QString fileBaseName() {
    return db_.databaseName().section('/', -1)
           .section('\\', -1);
  };

  bool beginTransaction();
  bool endTransaction();
  bool rollbackTransaction();

  QSqlQuery *newQuery() const;
  QSqlQuery *newQuery(const QString &sqlStr) const;
  QSqlQuery *newQuery(const QString &sqlStr, const QVariant &value) const;
  QSqlQuery *newQuery(const QString &sqlStr, const QVariantList &values) const;
  QSqlQuery *newTableQuery(const QString &table) const;
  QSqlQuery *newTableQuery(const QString &table, int id) const;
  QSqlQuery *newTableQuery(const QString &table, const QString &refField,
                           int id) const;
  QSqlQuery *newTableQuery(const QString &table, const QString &refField,
                           int id, const QString &sortField) const;

  QSqlQueryModel *newModel(const QString &sqlStr) const;
  QSqlQueryModel *newModel(const QString &sqlStr, const QVariant &value) const;
  QSqlQueryModel *newModel(const QString &sqlStr, 
                           const QVariantList &values) const;
  QSqlQueryModel *newTableModel(const QString &table) const;
  QSqlQueryModel *newTableModel(const QString &table, int id) const;
  QSqlQueryModel *newTableModel(const QString &table, const QString &refField,
                                int id) const;
  QSqlQueryModel *newTableModel(const QString &table, const QString &refField,
                                int id, const QString &sortField) const;
  QSqlQueryModel *newTableModel(const QString &table,
                                const QString &sortField) const;
  QSqlQueryModel *newTableModel(const QString &table,
                                const QMultiHash<QString, QStringList> &relations) const;

  void updateModel(QSqlQueryModel *model) const;

  inline int getModelId(QSqlQueryModel *model, int row) const { 
    return model->index(row, 0).data().toInt(); 
  }

  int getModelIdRow(QSqlQueryModel *model, int id) const;


  bool exist(const QString &table) const;
  bool exist(const QString &table, const QString &field,
             const QVariant &content) const;
  int existId(const QString &table, const QString &field,
              const QVariant &content) const;


  int getNumRows(const QString &table) const;
  int insertRow(const QString &table, const QHash<QString, QVariant> &values,
                bool ignore = false);
  bool updateRow(const QString &table, int id,
                 const QHash<QString, QVariant> &values) const;
  bool removeRow(const QString &table, int id);

  void emptyTables(const QList<QString> &tables) const;
  void emptyTablesNoKeys(const QList<QString> &tables) const;

  static const char *DatetimeFormat;

  inline QSqlError lastError() const {
    return db_.lastError();
  }

 private:
  Q_DISABLE_COPY(DB);

  int openImportDB(const QString &fileName, QSqlDatabase *db);
  void fetchAllData(QSqlQueryModel *model) const;

  QSqlDatabase db_;
  int nNestedTrans_;
};

} // namespace LoboLab
