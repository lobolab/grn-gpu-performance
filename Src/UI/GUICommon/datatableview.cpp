// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#include "datatableview.h"

#include <QHeaderView>
#include <QDebug>

namespace LoboLab {

DataTableView::DataTableView(QWidget * parent)
  : QTableView(parent) {
  sortProxyModel_ = new NaturalSortProxyModel(this);
  QTableView::setModel(sortProxyModel_);
}

DataTableView::~DataTableView() {
  delete sortProxyModel_->sourceModel();
}

void DataTableView::setModel(QAbstractItemModel *model, int sortColumn,
                             Qt::SortOrder order) {
  sortProxyModel_->setSourceModel(model);

  if (model) {
    fetchAllData();

    setColumnHidden(0, true); // 0 is Id column
    setSelectionBehavior(QAbstractItemView::SelectRows);
    setEditTriggers(QAbstractItemView::NoEditTriggers);
    setSortingEnabled(true);
    sortByColumn(sortColumn, order);
  }
}

void DataTableView::modelUpdated() {
  fetchAllData();
}

void DataTableView::fetchAllData() {
  //// We need to fetch all the data to order it properly
  QModelIndex invalidIndex;
  while (sortProxyModel_->canFetchMore(invalidIndex))
    sortProxyModel_->fetchMore(invalidIndex);
}

}