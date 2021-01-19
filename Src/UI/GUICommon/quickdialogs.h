// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#pragma once

#include <QWidget>

namespace LoboLab {

class DB;

namespace QuickDialogs {
  void getUniqueText(DB *db, QWidget *parent, const QString &title,
                     const QString &table, const QString &field, 
                     QString &newText);
} // namespace QuickDialogs

} // namespace LoboLab
