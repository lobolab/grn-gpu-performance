// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#include "mainwindow.h"

#include <QApplication>

#include <qwt_text.h>
#include <qwt_mathml_text_engine.h>

int main(int argc, char *argv[]) {
  //_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

  QApplication a(argc, argv);
  
  // This needs to be done only once before using the MathML engine
  QwtText::setTextEngine(QwtText::MathMLText, new QwtMathMLTextEngine());

  LoboLab::MainWindow mw;
  mw.show();

  return a.exec();
}
