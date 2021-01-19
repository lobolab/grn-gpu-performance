// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#pragma once

#include "Common/mathalgo.h"
#include "Model/model.h"

#include <QHash>
#include <QPoint>
#include <QRect>
#include <QVector>

namespace LoboLab {

class MorphologyImage;

class SimState {
  friend class SimCalculator;

 public:
  SimState(const QSize &size, int nInputMorphogens, int nTargetMorphogens, double distErrorThreshold);
  SimState(const SimState &source);
  SimState &operator=(const SimState &source);
  ~SimState();

  void initialize(int nProducts);
  void copyFrom(const SimState &source);

  inline Eigen::MatrixXd &product(int i) { return products_[i]; }
  inline const Eigen::MatrixXd &product(int i)  const { return products_.at(i);}
  inline double distErrorThreshold() const { return distErrorThreshold_; }
  inline int nProducts() const { return nProductsInUse_; }
  inline int nInputMorphogens() const { return nInputMorphogens_; }
  inline int nTargetMorphogens() const { return nTargetMorphogens_; }
  inline const QSize &size() const { return size_; }

  void clearProducts();

  void cloneCell(const QPoint &orig, const QPoint &dest);
  void deleteCell(int i, int j);

  double negativeLogLikelihood() const;
  double negativeLogLikelihood(const SimState &other) const;
  double calcDistSimple(const SimState &other) const;
  double calcDist(const SimState &other) const;

  void loadMorphologyImage(const MorphologyImage &morphologyImage, int nMorphogens);
  
 private:
  void adjustProductsUsed(int nProducts);

  double calcMinNearCellDist(int i1, int j1, int i2, int j2, 
                      const SimState &simState2) const;
  double calcProdDist2(int i1, int j1, int i2, int j2, 
                       const QList<Eigen::MatrixXd> &products2) const;

  QSize size_;
  int nInputMorphogens_;
  int nTargetMorphogens_;
  double distErrorThreshold_;
  int nProductsInUse_; // products can contain more matrices than necessary
  // Eigen Arrays are used transposed: width = rows, height = cols
  // and it is more efficient to traverse the Array by height and then by width,
  // since the default storage order is by column.
  QList<Eigen::MatrixXd> products_;
};

} // namespace LoboLab
