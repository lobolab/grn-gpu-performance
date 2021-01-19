// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#pragma once

#include <QList>

namespace LoboLab {

class Search;
class Individual;

class ErrorCalculator {

 public:
  ErrorCalculator();
  virtual ~ErrorCalculator();

  virtual void process(int nDeme, const QList<Individual*> &individuals) = 0;
  virtual int waitForAnyDeme() = 0;
};

} // namespace LoboLab
