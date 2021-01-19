// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#pragma once

#include <QDialog>
#include <QLineEdit>
#include <QLabel>
#include <QDialogButtonBox>
#include <QComboBox>
#include <QGroupBox>
#include <QSqlQueryModel>
#include <QHBoxLayout>

namespace LoboLab {

class Experiment;

class DB;

class ExperimentForm : public QDialog {
  Q_OBJECT

  public:
    ExperimentForm(Experiment *experiment, DB *db, QWidget *parent = NULL);
    virtual ~ExperimentForm();

    public slots:

    private slots :
      void nameChanged();
    void imageDeleted();

    void keyPressEvent(QKeyEvent *e);
    void formAccepted();

  private:
    void createWidgets();
    void readSettings();
    void writeSettings();

    QLineEdit *nameEdit_;
    QLineEdit *nameInputImage_;
    QLineEdit *nameOutputImage_;
    QString originalName_;

    Experiment *experiment_;
    DB *db_;

    QWidget *inputImageWidget; 
    QHBoxLayout *inputImageLayout_; 

    QWidget *outputImageWidget; 
    QHBoxLayout *outputImageLayout_; 

};

} // namespace LoboLab
