// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#pragma once

#include <QTextStream>

namespace LoboLab {

class ModelLink {
 public:
   // Random link
  ModelLink(int regulator = 0, int regulated = 0, bool isActivation = false); 
  ModelLink(int regulator, int regulated, double disConst, double hillCoef,
            bool isAndReg);
  virtual ~ModelLink();
  ModelLink(const ModelLink &source);

  inline int regulatorProdLabel() const { return regulator_; }
  inline int regulatedProdLabel() const { return regulated_; }
  inline double disConst() const { return disConst_; }
  inline double hillCoef() const { return hillCoef_; }
  inline bool isAndReg() const { return isAndReg_; }

  inline void setRegulator(int regulatorLabel) {regulator_ = regulatorLabel;}
  inline void setRegulated(int regulatedLabel) {regulated_ = regulatedLabel;}

  inline int complexity() const { return 1; }

  void mutateParams(int mutationProb);
  //void mutateParam();

  // Text Serialization
  void loadFromString(QString &str);
  QString toString();

  friend QTextStream &operator<<(QTextStream &stream, const ModelLink &prod);
  friend QTextStream &operator>>(QTextStream &stream, ModelLink &prod);

 private:
  int regulator_; // These are product labels
  int regulated_;
  double disConst_;
  double hillCoef_;
  bool isAndReg_; // if false, is OR regulation
};

} // namespace LoboLab
