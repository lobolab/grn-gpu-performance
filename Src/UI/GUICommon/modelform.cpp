// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#include "modelform.h"

#include "Model/model.h"
#include "Private/modelprodlistwidget.h"
#include "Private/modelformulawidget.h"
#include "modelgraphview.h"

#include <QVBoxLayout>
#include <QDialogButtonBox>
#include <QMessageBox>
#include <QSettings>
#include <QPushButton>
#include <QLabel>
#include <QInputDialog>
#include <QApplication>
#include <QClipboard>
#include <QSplitter>

namespace LoboLab {

ModelForm::ModelForm(Model *m, 
  DB *d , int nMorphogens, QWidget *parent, bool autoDelete,
  const QString &windowTitle)
  : QDialog(parent),
  model_(m),
  autoDelete_(autoDelete),
  nMorphogens_(nMorphogens),
  isUpdating_(false) {
  setWindowTitle(windowTitle);
  createWidgets(nMorphogens_);
  readSettings();
  updateLabels();
}

ModelForm::~ModelForm() {
  writeSettings();
  if (autoDelete_)
    delete model_;
}

void ModelForm::readSettings() {
  QSettings settings;
  settings.beginGroup("ModelForm");
  resize(settings.value("size", QSize(450, 600)).toSize());
  move(settings.value("pos", pos()).toPoint());
  if (settings.value("maximized", false).toBool())
    showMaximized();
}

void ModelForm::writeSettings() {
  QSettings settings;
  settings.beginGroup("ModelForm");
  if (isMaximized())
    settings.setValue("maximized", isMaximized());
  else {
    settings.setValue("maximized", false);
    settings.setValue("size", size());
    settings.setValue("pos", pos());
  }

  settings.endGroup();
}

void ModelForm::formAccepted() {
  accept();
}

Model *ModelForm::getModel() const {
  return model_;
}

void ModelForm::createWidgets(int nMorphogens) {
  setWindowFlags(Qt::Window);

  nProductsLabel_ = new QLabel();
  nLinksLabel_ = new QLabel();
  complexityLabel_ = new QLabel();

  QPushButton *addProductButton = new QPushButton("Add prod.");
  connect(addProductButton, SIGNAL(clicked()), this, SLOT(addProduct()));

  QPushButton *removeProductButton = new QPushButton("Remove prod.");
  connect(removeProductButton, SIGNAL(clicked()), this, SLOT(removeProduct()));


  modelProdListWidget_ = new ModelProdListWidget(model_, nMorphogens, this);
  connect(modelProdListWidget_, SIGNAL(modified()),
    this, SLOT(modelProdListWidgetChanged()));

  modelGraphView_ = new ModelGraphView(model_, nMorphogens);

  modelTextEdit_ = new QTextEdit();
  modelTextEdit_->setAcceptRichText(false);
  QString modelStr = model_->toString();
  modelTextEdit_->setText(modelStr);
  connect(modelTextEdit_, SIGNAL(textChanged()), this, SLOT(textChanged()));

  modelFormulaWidget_ = new ModelFormulaWidget(model_, nMorphogens);

  hideNotUsedCheckBox_ = new QCheckBox("Hide not used products");
  connect(hideNotUsedCheckBox_, SIGNAL(stateChanged(int)),
    this, SLOT(hideNotUsedCheckBoxChanged(int)));

  QPushButton *resetButton = new QPushButton("Clear");
  connect(resetButton, SIGNAL(clicked()), this, SLOT(clearModel()));

  QPushButton *copyMathMLButton = new QPushButton("Copy MathML");
  connect(copyMathMLButton, SIGNAL(clicked()), this, SLOT(copyMathML()));

  QDialogButtonBox *buttonBox = new QDialogButtonBox(QDialogButtonBox::Ok |
    QDialogButtonBox::Cancel);
  connect(buttonBox, SIGNAL(accepted()), this, SLOT(formAccepted()));
  connect(buttonBox, SIGNAL(rejected()), this, SLOT(reject()));

  QHBoxLayout *butLay = new QHBoxLayout();
  butLay->addWidget(resetButton);
  butLay->addWidget(copyMathMLButton);
  butLay->addStretch();
  butLay->addWidget(hideNotUsedCheckBox_);
  butLay->addWidget(buttonBox);


  QHBoxLayout *topLay = new QHBoxLayout();
  topLay->addWidget(new QLabel("Prod:"));
  topLay->addWidget(nProductsLabel_);
  topLay->addWidget(new QLabel("Links:"));
  topLay->addWidget(nLinksLabel_);
  topLay->addWidget(new QLabel("Complex:"));
  topLay->addWidget(complexityLabel_);
  topLay->addWidget(addProductButton);
  topLay->addWidget(removeProductButton);

  QSplitter *leftSplitter = new QSplitter(Qt::Vertical);
  leftSplitter->addWidget(modelFormulaWidget_);
  leftSplitter->addWidget(modelProdListWidget_);
  leftSplitter->addWidget(modelTextEdit_);
  leftSplitter->setSizes(QList<int>() << 10 << 500 << 10);

  QVBoxLayout *leftLay = new QVBoxLayout();
  leftLay->addLayout(topLay);
  leftLay->addWidget(leftSplitter, 1);

  QVBoxLayout *rightLay = new QVBoxLayout();
  rightLay->addWidget(modelGraphView_);
  rightLay->setMargin(0);

  QWidget *leftWidget = new QWidget();
  leftWidget->setLayout(leftLay);
  leftLay->setMargin(0);

  QWidget *rightWidget = new QWidget();
  rightWidget->setLayout(rightLay);

  QSplitter *centerSplitter = new QSplitter(Qt::Horizontal);
  centerSplitter->addWidget(leftWidget);
  centerSplitter->addWidget(rightWidget);
  centerSplitter->setSizes(QList<int>() << 100 << 400);

  QVBoxLayout *mainLay = new QVBoxLayout();
  mainLay->addWidget(centerSplitter, 1);
  mainLay->addLayout(butLay);

  setLayout(mainLay);
}

void ModelForm::updateLabels() {
  nProductsLabel_->setText(QString("%1 (%2)")
    .arg(model_->calcProductLabelsInUse(nMorphogens_).size())
    .arg(model_->nProducts()));
  nLinksLabel_->setText(QString("%1 (%2)")
    .arg(model_->calcLinksInUse(4).size())
    .arg(model_->nLinks()));
  complexityLabel_->setText(QString("%1 (%2)")
    .arg(model_->calcComplexityInUse(4))
    .arg(model_->calcComplexity()));
}

void ModelForm::modelProdListWidgetChanged() {
  isUpdating_ = true;
  modelTextEdit_->setText(model_->toString());
  modelGraphView_->updateModel(hideNotUsedCheckBox_->isChecked());
  updateLabels();
  isUpdating_ = false;
}

void ModelForm::textChanged() {
  if (!isUpdating_) {
    model_->loadFromString(modelTextEdit_->toPlainText());
    updateLabels();
    modelProdListWidget_->updateList(nMorphogens_);
    modelFormulaWidget_->updateFormula(nMorphogens_);
    modelGraphView_->updateModel(hideNotUsedCheckBox_->isChecked());
  }
}

void ModelForm::addProduct() {
}

void ModelForm::removeProduct() {
  bool ok;
  int delProd = QInputDialog::getInt(this, tr("Select Product to Remove"),
    tr("Product ind:"), 0, 0,
    model_->nProducts(), 1, &ok);

  if (ok) {
    model_->removeProduct(delProd);
    isUpdating_ = true;
    updateLabels();
    modelProdListWidget_->updateList(nMorphogens_);
    modelFormulaWidget_->updateFormula(nMorphogens_);
    modelGraphView_->updateModel(hideNotUsedCheckBox_->isChecked());
    modelTextEdit_->setText(model_->toString());
    isUpdating_ = false;
  }
}

void ModelForm::clearModel() {
  model_->clear();
  modelProdListWidget_->update();
  modelFormulaWidget_->updateFormula(nMorphogens_);
  modelTextEdit_->setText(model_->toString());
}

void ModelForm::hideNotUsedCheckBoxChanged(int state) {
  modelGraphView_->updateModel(state > 0);
}

void ModelForm::copyMathML() {
  QClipboard *clipboard = QApplication::clipboard();
  clipboard->setText(modelFormulaWidget_->mathMLStr());
}

}