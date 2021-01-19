// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#include "modelprodlistwidget.h"

#include "colorsconfig.h"
#include "Model/model.h"
#include "Model/modelprod.h"
#include "Model/modellink.h"
#include "Common/mathalgo.h"

#include <QAction>
#include <QHeaderView>
#include <QInputDialog>
#include <QMessageBox>
#include <QMenu>

namespace LoboLab {

ModelProdListWidget::ModelProdListWidget(Model *m, int nMorphogens, QWidget * parent)
  : QListWidget(parent), model_(m) {
  createActions();

  setContextMenuPolicy(Qt::CustomContextMenu);
  connect(this, SIGNAL(customContextMenuRequested(const QPoint &)),
    this, SLOT(contextMenuCalled(const QPoint &)));

  connect(this, SIGNAL(itemDoubleClicked(QListWidgetItem *, int)),
    this, SLOT(editItem(QListWidgetItem *)));

  createList(nMorphogens);
}

ModelProdListWidget::~ModelProdListWidget() {
}

void ModelProdListWidget::createActions() {
  // Context menu when clicking in a blank space
  blankContextMenu_ = new QMenu(this);

  QAction *action = new QAction("Add new product", blankContextMenu_);
  action->setIcon(QIcon(":/Images/graphicsRegion.png"));
  connect(action, SIGNAL(triggered()),
    this, SLOT(addModelProd()));
  blankContextMenu_->addAction(action);

  // Context menu when clicking in a product
  prodContextMenu_ = new QMenu(this);

  action = new QAction("Replace product", prodContextMenu_);
  action->setIcon(QIcon(":/Images/graphicsRegion.png"));
  connect(action, SIGNAL(triggered()),
    this, SLOT(replaceModelProd()));
  prodContextMenu_->addAction(action);

  action = new QAction("Remove component", prodContextMenu_);
  action->setIcon(QIcon(":/Images/graphicsRegion.png"));
  connect(action, SIGNAL(triggered()),
    this, SLOT(removeSelectedModelProd()));
  prodContextMenu_->addAction(action);
}

void ModelProdListWidget::contextMenuCalled(const QPoint &pos) {
  QListWidgetItem *item = itemAt(pos);
  if (item) {
    setCurrentItem(item);
    prodContextMenu_->popup(mapToGlobal(pos));
  }
  else
    blankContextMenu_->popup(mapToGlobal(pos));
}

void ModelProdListWidget::addModelProd() {
}

void ModelProdListWidget::replaceModelProd() {
}

ModelProd *ModelProdListWidget::createModelProd() {
  ModelProd *newProd = NULL;

  return newProd;
}

void ModelProdListWidget::removeSelectedModelProd() {
}

int ModelProdListWidget::selectProduct(bool *ok, const QString &label) {
  return QInputDialog::getInt(this, tr("Select Product"), label, 0, 0,
    model_->nProducts(), 1, ok);
}

double ModelProdListWidget::selectRealValue(bool *ok, const QString &label) {
  return QInputDialog::getDouble(this, tr("Select Real Value"), label, 0,
    -1000, 1000, 3, ok);
}

bool ModelProdListWidget::selectBoolValue(bool *ok, const QString &label) {
  return QInputDialog::getInt(this, tr("Select Yes or No"), label, 1,
    0, 1, 1, ok);
}

int ModelProdListWidget::selectDirection(bool *ok) {
  QStringList items;
  items << "North" << "South";
  QString selectedItem = QInputDialog::getItem(this, "Select direction",
    "Direction:", items, 0, false, ok);

  return items.indexOf(selectedItem);
}

void ModelProdListWidget::clearActions() {
  while (!actions().isEmpty())
    removeAction(actions().first());
}

void ModelProdListWidget::updateList(int nMorphogens) {
  clear();
  createList(nMorphogens);
}

void ModelProdListWidget::createList(int nMorphogens) {
  QList<int> labels = model_->calcProductLabels().toList();
  qSort(labels);

  int colorInd = 0;
  int nProducts = labels.size();
  // Part of GUI. did not implement to take Number of Morphogens from the database.
  QSet<int> labelsInUse = model_->calcProductLabelsInUse(nMorphogens);

  for (int i = 0; i < nProducts; ++i) {
    int label = labels.at(i);
    ModelProd *modelProd = model_->prodWithLabel(labels.at(i));

    QString preStr;
    if (labels.last() > 9 && label < 10)
      preStr = ' ';

    QString str = preStr + QString("[%1] ini=%2 lim=%3 deg=%4 dif=%5").arg(label)
      .arg(modelProd->init()).arg(modelProd->lim())
      .arg(modelProd->deg()).arg(modelProd->dif());

    QListWidgetItem *item = new QListWidgetItem(str);
    item->setData(Qt::UserRole, i);

    QColor frontColor;
    QColor backColor;

    if (labelsInUse.contains(label)) {
      frontColor = Qt::black;
      backColor = (ColorsConfig::colorsCells[colorInd++ % ColorsConfig::nColors]);
      backColor.setAlpha(255 - 0x80);
    }
    else {
      frontColor = Qt::gray;
      backColor = Qt::white;
    }

    QFont font("Monospace");
    font.setStyleHint(QFont::TypeWriter);

    item->setFont(font);
    item->setForeground(QBrush(frontColor));
    item->setBackground(QBrush(backColor));
    addItem(item);

    // Links
    QList<ModelLink*> links = model_->linksToLabel(label);

    if (labels.last() > 9)
      preStr = "    ";
    else
      preStr = "   ";

    int nLinks = links.size();
    for (int j = 0; j < nLinks; ++j) {
      ModelLink *link = links.at(j);
      int regulatorLabel = link->regulatorProdLabel();
      QString postStr;
      if (regulatorLabel < 10)
        postStr = ' ';

      str = preStr + QString("%1-- [%2] %3%4 dc=%5 hc=%6")
        .arg(link->hillCoef() >= 0 ? '<' : '|')
        .arg(regulatorLabel).arg(postStr)
        .arg(link->isAndReg() ? " and" : " or ")
        .arg(link->disConst())
        .arg(link->hillCoef());

      QListWidgetItem *item = new QListWidgetItem(str);
      item->setData(Qt::UserRole, j);
      item->setFont(font);
      item->setForeground(QBrush(frontColor));
      item->setBackground(QBrush(backColor));

      addItem(item);
    }
  }

  //sortItems();
}

}