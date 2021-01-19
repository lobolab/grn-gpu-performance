// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#include "quickdialogs.h"
#include "DB/db.h"

#include <QInputDialog>
#include <QMessageBox>

namespace LoboLab {

void QuickDialogs::getUniqueText(DB *db, QWidget *parent, const QString &title,
                                 const QString &table, const QString &field, 
                                 QString &newText) {
  QString tableLow = table.toLower();
  QString fieldLow = field.toLower();
  bool ok = true;
  bool nameUnique = false;
  while (ok && !nameUnique) {
    QString labelLine = table + " " + fieldLow + ":";
    newText = QInputDialog::getText(parent, title, labelLine,
                                    QLineEdit::Normal, newText, &ok);

    if (ok) {
      if (newText.isEmpty())
        ok = false;
      else if (db->exist(table, field, newText)) {
        QMessageBox::warning(parent, table + " with same " +
                             fieldLow + " found", "The database already contains a " +
                             tableLow + " with the same " + fieldLow + ". "
                             "Please, specify a different " + fieldLow + " for this " +
                             tableLow + ".", QMessageBox::Ok, QMessageBox::Ok);
      } else
        nameUnique = true;
    }
  }

  if (!ok)
    newText.clear();
}


}