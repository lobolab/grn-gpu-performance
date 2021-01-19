// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#pragma once

#include <QDialog>
#include <QDialogButtonBox>
#include <QTextEdit>
#include <QLabel>
#include <QCheckBox>

#include <qwt_text_label.h>

namespace LoboLab {

class Model;
class ModelProdListWidget;
class ModelFormulaWidget;
class ModelGraphView;
class Search;
class DB;

class ModelForm : public QDialog {
  Q_OBJECT

  public:
    ModelForm(Model *model, DB *db, int nMorphogens, QWidget *parent = NULL, bool autoDelete = false,
      const QString &windowTitle = "Edit model");
    virtual ~ModelForm();

    Model *getModel() const;

    private slots:
    void formAccepted();
    void modelProdListWidgetChanged();
    void textChanged();
    void addProduct();
    void removeProduct();
    void clearModel();
    void hideNotUsedCheckBoxChanged(int state);
    void copyMathML();

  private:
    void readSettings();
    void writeSettings();
    void createWidgets(int nMorphogens);
    void updateLabels();

    Model *model_;

    QLabel *nProductsLabel_;
    QLabel *nLinksLabel_;
    QLabel *complexityLabel_;
    ModelProdListWidget *modelProdListWidget_;
    ModelFormulaWidget *modelFormulaWidget_;
    ModelGraphView *modelGraphView_;
    QTextEdit *modelTextEdit_;
    QCheckBox *hideNotUsedCheckBox_;

    bool autoDelete_;
    bool isUpdating_;

    int nMorphogens_;
};

} // namespace LoboLab
