// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#include "modelprod.h"
#include "model.h"
#include "Common/mathalgo.h"

namespace LoboLab {

ModelProd::ModelProd() 
  : label_(-1) {
}

ModelProd::ModelProd(int nMorphogens, int label)
    : label_(label) {
  init_ = zeroInit(nMorphogens);
  lim_ = MathAlgo::rand01();
  deg_ = MathAlgo::rand01();

  if (MathAlgo::rand100() < 10)
    dif_ = MathAlgo::rand01();
  else
    dif_ = 0;
}

ModelProd::ModelProd(int label, double init, double lim, double deg, double dif)
  : label_(label), init_(init), lim_(lim), deg_(deg), dif_(dif) {
}


ModelProd::~ModelProd() {
}

ModelProd::ModelProd(const ModelProd &source)
  : label_(source.label_), init_(source.init_), lim_(source.lim_), 
    deg_(source.deg_), dif_(source.dif_) {
}

int ModelProd::complexity() const {
  return 1 + (dif_ > 0);
}

void ModelProd::mutateParams(int mutationProb, int nMorphogens) {
  // Change lim
  if (MathAlgo::rand100() < mutationProb) 
    lim_ = MathAlgo::rand01();

  // Change deg
  if (MathAlgo::rand100() < mutationProb) 
    deg_ = MathAlgo::rand01();

  // Change diffusion
  if (MathAlgo::rand100() < mutationProb) { 
    if (MathAlgo::randBool())
      dif_ = MathAlgo::rand01();
    else
      dif_ = 0;
  }
}

double ModelProd::zeroInit(int nMorphogenes) const {
	if (label_ == 0 && nMorphogenes > 0) // red
		return -1;
	else if (label_ == 1 && nMorphogenes > 1) // green
		return -2;
	else if (label_ == 2 && nMorphogenes > 2) // blue
		return -3;
	else 
		return 0;
}


// Serialization
QTextStream &operator<<(QTextStream &stream, const ModelProd &prod) {
  stream << prod.label_ << ' ' << prod.init_ << ' ' << prod.lim_ << ' ' << 
    prod.deg_ << ' ' << prod.dif_;

  return stream;
}

QTextStream &operator>>(QTextStream &stream, ModelProd &prod) {
  stream >> prod.label_;

  char c;
  stream >> c;
  Q_ASSERT(c == ' ');
  
  prod.init_ = Model::parseDouble(stream);

  stream >> c;
  Q_ASSERT(c == ' ');

  prod.lim_ = Model::parseDouble(stream);

  stream >> c;
  Q_ASSERT(c == ' ');

  prod.deg_ = Model::parseDouble(stream);

  stream >> c;
  Q_ASSERT(c == ' ');

  prod.dif_ = Model::parseDouble(stream);

  return stream;
}


bool prodLabelLessThan(const ModelProd *p1, const ModelProd *p2) {
  return p1->label() < p2->label();
}

}