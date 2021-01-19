// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#pragma once

#include <QTextStream>

namespace LoboLab {

class ModelProd {
 public:
  explicit ModelProd(); // Empty product
  explicit ModelProd(int nMorphogens, int label); // Random product
  ModelProd(int label, double init, double lim, double deg, double dif);
  virtual ~ModelProd();
  ModelProd(const ModelProd &source);

  inline int label() const { return label_; }
  inline double init() const { return init_; }
  inline double lim() const { return lim_; }
  inline double deg() const { return deg_; }
  inline double dif() const { return dif_; }
  inline void setLim(double lim) { lim_ = lim; }
  inline void setInit(double init) { init_ = init; }
  
  inline void setLabel(int label) { label_ = label; }

  int complexity() const;

  void mutateParams(int mutationProb, int nMorphogens);

  // Text Serialization
  void loadFromString(QString &str);
  QString toString();

  friend QTextStream &operator<<(QTextStream &stream, const ModelProd &prod);
  friend QTextStream &operator>>(QTextStream &stream, ModelProd &prod);

 private:
  //double randInit(int nMorphogenes) const;
  double zeroInit(int nMorphogenes) const;

  int label_;
  double init_;
  double lim_;
  double deg_;
  double dif_;
};

bool prodLabelLessThan(const ModelProd *p1, const ModelProd *p2);

} // namespace LoboLab
