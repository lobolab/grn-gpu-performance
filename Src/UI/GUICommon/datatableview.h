// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#pragma once

#include "Private/naturalstringcompare.h"

#include <QTableView>
#include <QSortFilterProxyModel>
#include <QDebug>

namespace LoboLab {

class DataTableView : public QTableView {
 Q_OBJECT

 public:
  DataTableView(QWidget * parent = NULL);
  virtual ~DataTableView();

  virtual void setModel(QAbstractItemModel * model, int sortColumn = 1,
                        Qt::SortOrder order = Qt::AscendingOrder);

  inline QModelIndex currentIndex() const {
    return sortProxyModel_->mapToSource(QTableView::currentIndex());
  }

  inline QModelIndex mapToSource(const QModelIndex &index) const {
    return sortProxyModel_->mapToSource(index);
  }

  void modelUpdated();

 public slots:
  inline void	setCurrentIndex(const QModelIndex &index) {
    QTableView::setCurrentIndex(sortProxyModel_->mapFromSource(index));
  }

  inline void	selectRow(int row) {
    QTableView::setCurrentIndex(sortProxyModel_->mapFromSource(
                                  sortProxyModel_->sourceModel()->index(row, 1)));
  }

  inline void	selectFirstRow() {
    QTableView::setCurrentIndex(sortProxyModel_->index(0, 0));
  }


 protected:
  void dataChanged(const QModelIndex & topLeft, 
                   const QModelIndex & bottomRight) {
    QTableView::dataChanged(topLeft, bottomRight);
    qDebug() << "dataChanged";
  }

 private:
  class NaturalSortProxyModel : public QSortFilterProxyModel {
   public:
    NaturalSortProxyModel(QObject * parent = 0) 
        : QSortFilterProxyModel(parent) {
      setDynamicSortFilter(true);
    }

    QVariant headerData(int section, Qt::Orientation orientation,
                        int role) const {
      if (role != Qt::DisplayRole)
        return QVariant();
      if (orientation == Qt::Horizontal)
        return QSortFilterProxyModel::headerData(section, orientation, role);
      return section + 1;
    }

   protected:
    bool lessThan(const QModelIndex &left, const QModelIndex &right) const {
      return NaturalStringCompare::naturalStringCaseInsensitiveCompareLessThan(
               sourceModel()->data(left).toString(),
               sourceModel()->data(right).toString());
    }
  };

  void fetchAllData();

  NaturalSortProxyModel *sortProxyModel_;
};

} // namespace LoboLab
