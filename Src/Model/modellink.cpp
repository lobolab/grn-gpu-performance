// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#include "modellink.h"
#include "model.h"
#include "Common/mathalgo.h"

namespace LoboLab {

ModelLink::ModelLink(int rtor, int rted, bool isActivation)
  : regulator_(rtor), regulated_(rted) {
  disConst_ = MathAlgo::rand01();
  if (isActivation)
    hillCoef_ = MathAlgo::rand15();
  else
    if (MathAlgo::randBool())
      hillCoef_ = MathAlgo::rand15();
    else
      hillCoef_ = -1.0 * MathAlgo::rand15();

  isAndReg_ = (hillCoef_ < 0) || MathAlgo::randBool();
}

ModelLink::ModelLink(int rtor, int rted, double d, double h, bool i)
  : regulator_(rtor), regulated_(rted), disConst_(d), hillCoef_(h), 
    isAndReg_(i) {
}


ModelLink::~ModelLink() {}

ModelLink::ModelLink(const ModelLink &source)
  : regulator_(source.regulator_), regulated_(source.regulated_),
    disConst_(source.disConst_), hillCoef_(source.hillCoef_),
    isAndReg_(source.isAndReg_) {
}

void ModelLink::mutateParams(int mutationProb) {
  if (MathAlgo::rand100() < mutationProb)
    disConst_ = MathAlgo::rand01();

  if (MathAlgo::rand100() < mutationProb) {
    if (MathAlgo::randBool())
      hillCoef_ = MathAlgo::rand15();
    else
      hillCoef_ = -1.0 * MathAlgo::rand15();
  }

  if (hillCoef_ < 0)
    isAndReg_ = true;
  else if (MathAlgo::rand100() < mutationProb)
    isAndReg_ = MathAlgo::randBool();
}

// Serialization
QTextStream &operator<<(QTextStream &stream, const ModelLink &prod) {
  stream << prod.regulator_ << ' ' << prod.regulated_ << ' ' << prod.disConst_ <<
         ' ' << prod.hillCoef_ << ' ' << prod.isAndReg_;

  return stream;
}

QTextStream &operator>>(QTextStream &stream, ModelLink &prod) {
  stream >> prod.regulator_;

  char c;
  stream >> c;
  Q_ASSERT(c == ' ');

  stream >> prod.regulated_;

  stream >> c;
  Q_ASSERT(c == ' ');

  prod.disConst_ = Model::parseDouble(stream);

  stream >> c;
  Q_ASSERT(c == ' ');

  prod.hillCoef_ = Model::parseDouble(stream);

  stream >> c;
  Q_ASSERT(c == ' ');

  int b;
  stream >> b;
  prod.isAndReg_ = b;

  return stream;
}


}