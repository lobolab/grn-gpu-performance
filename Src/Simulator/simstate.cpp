// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#include "simstate.h"

#include "simulatorconfig.h"
#include "Common/log.h"
#include "Model/modelprod.h"
// #include "Experiment/manipulation.h"
// #include "Experiment/morphology.h"
#include "Experiment/morphologyimage.h"

#include <QImage>

namespace LoboLab {

SimState::SimState(const QSize &s, int nInputMorphogens, int nTargetMorphogens, double distErrorThreshold)
    : size_(s), nInputMorphogens_(nInputMorphogens), nTargetMorphogens_(nTargetMorphogens), distErrorThreshold_(distErrorThreshold), nProductsInUse_(0) {
  initialize(nInputMorphogens_ + nTargetMorphogens_);
}

SimState::~SimState() {
}

SimState::SimState(const SimState &source)
    : size_(source.size_), 
      nInputMorphogens_(source.nInputMorphogens_),
      nTargetMorphogens_(source.nTargetMorphogens_) {
  nProductsInUse_ = source.nProductsInUse_;
  for (int i = 0; i < nProductsInUse_; ++i)
    products_.append(source.products_[i]);
}

SimState &SimState::operator=(const SimState &source) {
  Q_ASSERT(size_ == source.size_);
  nInputMorphogens_ = source.nInputMorphogens_;
  nTargetMorphogens_ = source.nTargetMorphogens_;

  adjustProductsUsed(source.nProductsInUse_);
  for (int i = 0; i < nProductsInUse_; ++i)
    products_[i] = source.products_[i];

  return *this;
}

void SimState::copyFrom(const SimState &source) {
  Q_ASSERT(size_ == source.size_);

  if (nProductsInUse_ < source.nProductsInUse_)
    adjustProductsUsed(source.nProductsInUse_);
  
  for (int i = 0; i < source.nProductsInUse_; ++i)
    products_[i] = source.products_[i];
}

// Add matrices if there is less than nProducts
void SimState::adjustProductsUsed(int nProducts) {
  nProductsInUse_ = nProducts;
  for (int i = products_.size(); i < nProductsInUse_; ++i)
    products_.append(Eigen::MatrixXd(size_.width(), size_.height()));
}

// Adjust and initialize to zero new products
void SimState::initialize(int nProducts) {
  int oldProducts = nProductsInUse_;
  adjustProductsUsed(nProducts);
  for (int i = oldProducts; i < nProductsInUse_; ++i)
    products_[i].setZero();
}

void SimState::loadMorphologyImage(const MorphologyImage &morphologyImage, int nMorphogens) {
  // Note: the image is upside-down, because in Qt the y-axis grows down, the opposite than in a Morphology
  // The Y-axis is reversed in a QImage with respect a Morphology
  if (nMorphogens > nProductsInUse_)
    adjustProductsUsed(nMorphogens);

  const QImage &image = morphologyImage.image();
  int imgWidth = image.width();
  int imgHeight = image.height();
  for (int j = 0; j<imgHeight; ++j) {
    QRgb* scanLine = (QRgb*)image.scanLine(j);
    
    for (int i = 0; i<imgWidth; ++i) {
      QRgb color = scanLine[i];
      for (int k = 0; k < nMorphogens; ++k)
        product(k)(i, imgHeight - j - 1) =
          ((double)((color >> (16 - (8 * k))) & 0xff)) / 255.0;
    }
  }
}

void SimState::clearProducts() {
  for (int k=0; k<nProductsInUse_; ++k)
    products_[k].setZero();
}

void SimState::cloneCell(const QPoint &orig, const QPoint &dest) {
  for (int k=0; k<nProductsInUse_; ++k) {
    products_[k](dest.x(), dest.y()) =
      products_[k](orig.x(), orig.y());
  }
}

void SimState::deleteCell(int i, int j) {
    for (int k=0; k<nProductsInUse_; ++k) {
      products_[k](i, j) = 0;
    }
}

// Calculate Log Likelihood between this state and zero for nProductsInUse.
double SimState::negativeLogLikelihood() const {
  const double noCellDif = nProductsInUse_;//sqrt((double)nInputMorphogens);
  double dist;

  dist = 0;
  for (int i = 0; i < size_.width(); ++i) {
    for (int j = 0; j < size_.height(); ++j) {
      for (int k = 0; k < nProductsInUse_; ++k) {
        double a = products_[k](i, j);

        double sub = a - distErrorThreshold_;
        if (sub > 0)
          dist += log(1 + sub);
      }
    }
  }

  // Distance is divided by number of cells to normalize.
  dist /= size_.width() * size_.height();

  return dist;
}

// Calculate Log Likelihood between two states for nTargetMorphogens
double SimState::negativeLogLikelihood(const SimState &simState2) const {
	double dist = 0;
	const QList<Eigen::MatrixXd> products2 = simState2.products_;

	for (int i = 0; i < size_.width(); ++i) {
		for (int j = 0; j < size_.height(); ++j) {
			for (int k = 0; k < nTargetMorphogens_; ++k) {
        double a = products_[k + nInputMorphogens_](i, j);
        double b = products2[k](i, j);

        double absSub = fabs(a - b) - distErrorThreshold_;
        if (absSub > 0) {
          double d = log(1 + absSub);
          dist += d;
        }
			}
	  }
	}

	// Distance is divided by domain size to normalize.
	dist /= size_.width() * size_.height();
  
	return dist;
}

// Calculate Log Likelihood between two states for nTargetMorphogens
double SimState::calcDistSimple(const SimState &simState2) const {
  const double noCellDif = nTargetMorphogens_;//sqrt((double)nInputMorphogens);
  double dist;

  const QList<Eigen::MatrixXd> products2 = simState2.products_;

  dist = 0;
  double totalConc = 0;
  for (int i = 0; i < size_.width(); ++i) {
    for (int j = 0; j < size_.height(); ++j) {
      for (int k = 0; k < nTargetMorphogens_; ++k) {
        if (products_[k + nInputMorphogens_](i, j) > 0.9)
          totalConc += products_[k + nInputMorphogens_](i, j);
      }
    }
    double absSub = fabs(totalConc - 1000) - distErrorThreshold_;
    if (absSub > 0)
      dist += log(1 + absSub);
  }
  // Distance is divided by number of cells to normalize.
  dist /= size_.width() * size_.height();

  return dist;
}

// Distance between two sim states
double SimState::calcDist(const SimState &simState2) const {
  const double noCellDif = nInputMorphogens_;
  double dist;
  dist = 0;
  for (int i = 0; i < size_.width(); ++i) {
    for (int j = 0; j < size_.height(); ++j) {
      dist += calcMinNearCellDist(i, j, i, j, simState2);
    }
  }

  dist /= size_.width() * size_.height();

  return dist;
}

double SimState::calcMinNearCellDist(int i1, int j1, int i2, int j2, 
                                     const SimState &simState2) const {
  double minProdDist2 = calcProdDist2(i1, j1, i2, j2, simState2.products_);
  if (minProdDist2 == 0)
    return 0;

  if (j2+1 < size_.height()) {
    double d = calcProdDist2(i1, j1, i2, j2+1, simState2.products_);
    if (d < minProdDist2) {
      minProdDist2 = d;
      
      if (minProdDist2 == 0)
        return 0;
    }
  }
    
  if (j2-1 >= 0) {
    double d = calcProdDist2(i1, j1, i2, j2-1, simState2.products_);
    if (d < minProdDist2) {
      minProdDist2 = d;
      
      if (minProdDist2 == 0)
        return 0;
    }
  }
    
  if (i2+1 < size_.width()) {
    double d = calcProdDist2(i1, j1, i2+1, j2, simState2.products_);
    if (d < minProdDist2) {
      minProdDist2 = d;
      
      if (minProdDist2 == 0)
        return 0;
    }
  }

  if (i2-1 >= 0) {
    double d = calcProdDist2(i1, j1, i2-1, j2, simState2.products_);
    if (d < minProdDist2) {
      minProdDist2 = d;
      
      if (minProdDist2 == 0)
        return 0;
    }
  }

  




  if (j2+2 < size_.height()) {
    double d = calcProdDist2(i1, j1, i2, j2+2, simState2.products_);
    if (d < minProdDist2) {
      minProdDist2 = d;
      
      if (minProdDist2 == 0)
        return 0;
    }
  }
    
  if (j2-2 >= 0) {
    double d = calcProdDist2(i1, j1, i2, j2-2, simState2.products_);
    if (d < minProdDist2) {
      minProdDist2 = d;
      
      if (minProdDist2 == 0)
        return 0;
    }
  }
    
  if (i2+2 < size_.width()) {
    double d = calcProdDist2(i1, j1, i2+2, j2, simState2.products_);
    if (d < minProdDist2) {
      minProdDist2 = d;
      
      if (minProdDist2 == 0)
        return 0;
    }
  }

  if (i2-2 >= 0) {
    double d = calcProdDist2(i1, j1, i2-2, j2, simState2.products_);
    if (d < minProdDist2) {
      minProdDist2 = d;
      
      if (minProdDist2 == 0)
        return 0;
    }
  }

  // Log dist
  return log(1 + sqrt(minProdDist2));
}

// The result is squared
double SimState::calcProdDist2(int i1, int j1, int i2, int j2, 
                               const QList<Eigen::MatrixXd> &products2) const {
  double prodDist = 0;
  for (int k = 0; k < nInputMorphogens_; ++k) {
    double a = products_[k](i1, j1);
    double b = products2[k](i2, j2);

    // Threshold:
    double absSub = fabs(a - b) - distErrorThreshold_;
    
    if (absSub > 0)
      prodDist += pow(absSub, 2);
  }

  return prodDist;
}

}