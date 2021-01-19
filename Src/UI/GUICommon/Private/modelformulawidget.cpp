// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#include "modelformulawidget.h"

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

ModelFormulaWidget::ModelFormulaWidget(Model *m, int nMorphogens, QWidget * parent)
  : FormulaWidget(parent),
  model_(m) {
  updateFormula(nMorphogens);
}

ModelFormulaWidget::~ModelFormulaWidget() {
}

void ModelFormulaWidget::updateFormula(int nMorphogens) {
  QString mathMLStr("<math xmlns='http://www.w3.org/1998/Math/MathML'>");
  QList<int> usedLabels = model_->calcProductLabelsInUse(nMorphogens).toList();
  qSort(usedLabels);
  int nProducts = usedLabels.size();

  QStringList names;

  int n = nProducts - names.size();
  for (int i = 0; i < n; ++i)
    names.append(QChar(i + 97));

  mathMLStr += "<mtable>";

  for (int i = 0; i < nProducts; ++i) {
    ModelProd *prod = model_->prodWithLabel(usedLabels[i]);

    mathMLStr += "<mtr><mtd columnalign='right'>";
    mathMLStr += QString("<mfrac><mrow><mo>&#x2202;</mo><mi mathvariant='italic'>%1</mi></mrow><mrow><mo>&#x2202;</mo><mi>t</mi></mrow></mfrac>").arg(names.at(i));
    mathMLStr += "</mtd><mtd columnalign='center'><mo>=</mo></mtd><mtd columnalign='left'>";
    
    // links
    QList<ModelLink*> links = model_->linksToLabel(prod->label());
    qSort(links.begin(), links.end(), linkRegulatorLessThan);
    int nLinks = links.size();

    if (nLinks == 1 && links[0]->hillCoef() >= 0) {
      mathMLStr += QString("<mn>%1</mn>").arg(prod->lim(), 0, 'g', 2);
      mathMLStr += QString("<mo>&#x22c5;</mo>");

      ModelLink *link = links[0];
      mathMLStr += createHillFormula(link,
        names[usedLabels.indexOf(link->regulatorProdLabel())]);
    }
    else {
      QList<ModelLink*> orLinks, andLinks;
      bool anyActivator = false;
      for (int j = 0; j < nLinks; ++j) {
        ModelLink *link = links[j];
        if (link->isAndReg())
          andLinks.append(link);
        else
          orLinks.append(link);

        if (link->hillCoef() >= 0)
          anyActivator = true;
      }

      if (anyActivator) {
        mathMLStr += QString("<mn>%1</mn>").arg(prod->lim(), 0, 'g', 2);
        mathMLStr += QString("<mo>&#x22c5;</mo>");

        int nAndLinks = andLinks.size();
        int nOrLinks = orLinks.size();

        if (nAndLinks > 0) {
          mathMLStr += "<mi>min</mi><mfenced><mrow>"; // <mrow> because bug in word with commas in the middle with mfenced

          for (int j = 0; j < nAndLinks; ++j) {
            mathMLStr += createHillFormula(andLinks[j],
              names[usedLabels.indexOf(andLinks[j]->regulatorProdLabel())]);
            mathMLStr += "<mo>,</mo>";
          }
        }

        if (nOrLinks > 1)
          mathMLStr += "<mrow><mi>max</mi><mfenced><mrow>";

        for (int j = 0; j < nOrLinks; ++j) {
          mathMLStr += createHillFormula(orLinks[j],
            names[usedLabels.indexOf(orLinks[j]->regulatorProdLabel())]);
          if (j < nOrLinks - 1)
            mathMLStr += "<mo>,</mo>";
        }

        if (nOrLinks > 1)
          mathMLStr += "</mrow></mfenced></mrow>";

        if (nAndLinks > 0)
          mathMLStr += "</mrow></mfenced>";
      }
    }

    mathMLStr += QString("<mo>-</mo><mn>%1</mn><mo>&#x22c5;</mo><mi mathvariant='italic'>%2</mi>").arg(prod->deg(), 0, 'g', 2).arg(names.at(i));

    if (prod->dif() > 0)
      mathMLStr += QString("<mo>+</mo><mn>%1</mn><mo>&#x22c5;</mo><msup><mo>&#x2207;</mo><mn>2</mn></msup><mi mathvariant='italic'>%2</mi>").arg(prod->dif(), 0, 'g', 2).arg(names.at(i));

    mathMLStr += "</mtd></mtr>";
  }

  mathMLStr += "</mtable>";
  mathMLStr += "</math>";
  setMathMLStr(mathMLStr);
}

QString ModelFormulaWidget::createHillFormula(const ModelLink *link,
  const QString &regulator) const {
  double exp = abs(link->hillCoef());
  QString regTerm = QString("<msup><mi mathvariant='italic'>%1</mi>"
    "<mn>%2</mn></msup>").arg(regulator).arg(exp, 0, 'g', 2);
  QString constTerm = QString("<msup><mn>%1</mn><mn>%2</mn></msup>")
    .arg(link->disConst(), 0, 'g', 2).arg(exp, 0, 'g', 2);

  QString formulaStr;

  formulaStr += "<mfrac>";

  if (link->hillCoef() > 0)
    formulaStr += regTerm;
  else
    formulaStr += constTerm;

  formulaStr += "<mrow>" + constTerm + "<mo>+</mo>" + regTerm + "</mrow>";

  formulaStr += "</mfrac>";

  return formulaStr;
}

// For sorting
bool ModelFormulaWidget::linkRegulatorLessThan(ModelLink *l1, ModelLink *l2) {
  return (l1->regulatorProdLabel() < l2->regulatorProdLabel());
}

}