// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#include "experimentform.h"
#include "Experiment/experiment.h"
#include "Experiment/morphologyimage.h" 
#include "Private/morphologyimagewidget.h" 
#include "DB/db.h"
#include <QDataWidgetMapper>
#include <QGridLayout>
#include <QPushButton>
#include <QInputDialog>
#include <QKeyEvent>
#include <QMessageBox>
#include <QSettings>
#include <QSplitter>
#include <QFormLayout>

namespace LoboLab {

ExperimentForm::ExperimentForm(Experiment *e, DB *d, QWidget *parent)
  : QDialog(parent), db_(d), experiment_(e) {
  if (experiment_->id()) {
    setWindowTitle(tr("Edit experiment"));
    originalName_ = experiment_->name();
  }
  else
    setWindowTitle(tr("New experiment"));

  createWidgets();
  readSettings();
}

ExperimentForm::~ExperimentForm() {
  writeSettings();
}

void ExperimentForm::readSettings() {
  QSettings settings;
  settings.beginGroup("ExperimentForm");
  resize(800, 500);
  if (settings.value("maximized", false).toBool())
    showMaximized();
}

void ExperimentForm::writeSettings() {
  QSettings settings;
  settings.beginGroup("ExperimentForm");
  if (isMaximized())
    settings.setValue("maximized", isMaximized());
  else {
    settings.setValue("maximized", false);
    settings.setValue("size", size());
  }
  settings.endGroup();
}

void ExperimentForm::formAccepted() {
  QString name = experiment_->name();

  if ((!experiment_->id() || name != originalName_)
    && db_->exist("Experiment", "Name", name)) {
    QMessageBox::warning(this, "Experiment with same name found", "The database "
      "already contains an experiment with the same name. Please, "
      "specify a different name for this experiment.", QMessageBox::Ok,
      QMessageBox::Ok);
    nameEdit_->setFocus();
  }
  else 
    accept();
}

void ExperimentForm::createWidgets() {
 
  setWindowFlags(Qt::Window);


  nameEdit_ = new QLineEdit();
  nameEdit_->setText(experiment_->name());
  QLabel *nameLabel = new QLabel("Name:");
  nameLabel->setBuddy(nameEdit_);

  // Input Image
  nameInputImage_ = new QLineEdit();
  nameInputImage_->setText(experiment_->inputMorphology()->name());
  QLabel *nameLabelInputImage = new QLabel("Name:");
  nameLabelInputImage->setBuddy(nameInputImage_);

  inputImageLayout_ = new QHBoxLayout();
  inputImageWidget = new QWidget(this);
  if (experiment_->inputMorphology()) {
    MorphologyImageWidget *inputImgWidget =
      new MorphologyImageWidget(experiment_->inputMorphology(), this);

    connect(inputImgWidget, SIGNAL(deleted()), this, SLOT(imageDeleted()));
    inputImageLayout_->addWidget(inputImgWidget);
    inputImageLayout_->setMargin(0);

    inputImageWidget->setLayout(inputImageLayout_);
  }
  else
    inputImageWidget->setVisible(false);

  QGroupBox *inputImageBox = new QGroupBox("Input Image");

  // Output Image
  nameOutputImage_ = new QLineEdit();
  nameOutputImage_->setText(experiment_->outputMorphology()->name());
  connect(nameOutputImage_, SIGNAL(editingFinished()), this, SLOT(nameChanged()));
  QLabel *nameLabelOutputImage = new QLabel("Name:");
  nameLabelOutputImage->setBuddy(nameOutputImage_);

  outputImageLayout_ = new QHBoxLayout();
  outputImageWidget = new QWidget(this);
  if (experiment_->outputMorphology()) {
    MorphologyImageWidget *outputImgWidget =
      new MorphologyImageWidget(experiment_->outputMorphology(), this);

    connect(outputImgWidget, SIGNAL(deleted()), this, SLOT(imageDeleted()));
    outputImageLayout_->addWidget(outputImgWidget);
    outputImageLayout_->setMargin(0);

    outputImageWidget->setLayout(outputImageLayout_);
  }
  else
    outputImageWidget->setVisible(false);

  QGroupBox *outputImageBox = new QGroupBox("Output Image");

  QDialogButtonBox *buttonBox = new QDialogButtonBox(QDialogButtonBox::Ok |
    QDialogButtonBox::Cancel);
  connect(buttonBox, SIGNAL(accepted()), this, SLOT(formAccepted()));
  connect(buttonBox, SIGNAL(rejected()), this, SLOT(reject()));

  // Layouts

  QHBoxLayout *nameLay = new QHBoxLayout();
  nameLay->addWidget(nameLabel);
  nameLay->addWidget(nameEdit_);

  QGridLayout *inputImageLayout = new QGridLayout;
  inputImageLayout->addWidget(nameLabelInputImage, 0, 0, 1, 1);
  inputImageLayout->addWidget(nameInputImage_, 0, 1, 1, 1);
  inputImageLayout->addWidget(inputImageWidget, 1, 0, 1, 2);
  inputImageBox->setLayout(inputImageLayout);


  QGridLayout *outputImageLayout = new QGridLayout;
  outputImageLayout->addWidget(nameLabelOutputImage, 0, 0);
  outputImageLayout->addWidget(nameOutputImage_, 0, 1);
  outputImageLayout->addWidget(outputImageWidget, 1, 0, 1, 2);
  outputImageBox->setLayout(outputImageLayout);

  QHBoxLayout *imagelayout = new QHBoxLayout();
  imagelayout->addWidget(inputImageBox);
  imagelayout->addWidget(outputImageBox);

  QVBoxLayout *mainLay = new QVBoxLayout();
  mainLay->addLayout(nameLay);
  mainLay->addSpacing(10);
  mainLay->addLayout(imagelayout);
  mainLay->addSpacing(10);
  mainLay->addWidget(buttonBox);
  setMinimumSize(800, 500);

  setLayout(mainLay);
}

// private slot
void ExperimentForm::nameChanged() {
  experiment_->setName(nameEdit_->text());
}

// private slot
void ExperimentForm::imageDeleted() {
  if (!experiment_->inputMorphology())
    inputImageWidget->setVisible(false);
  if (!experiment_->outputMorphology())
    outputImageWidget->setVisible(false);
}

void ExperimentForm::keyPressEvent(QKeyEvent *e) {
  if (e->key() != Qt::Key_Escape)
    QDialog::keyPressEvent(e);
}
}
