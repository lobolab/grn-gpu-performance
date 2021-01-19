// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#include "concentimagethread.h"
#include "simulatorthread.h"
#include "Simulator/simulator.h"
#include "Simulator/simstate.h"
#include "../Private/colorsconfig.h"
#include <QImage>

namespace LoboLab {

ConcentImageThread::ConcentImageThread(Simulator &sim,
                                       double concentImageFactor,
                                       SimulatorThread *parent)
  : QThread(parent), simulator_(sim), size_(sim.domainSize()),
    image_(size_, QImage::Format_RGB32), factor_(concentImageFactor) {
  calcImage();
}


ConcentImageThread::~ConcentImageThread() {
  stopThread();
}

void ConcentImageThread::run() {
  mutex_.lock();
  endThread_ = false;

  while (!endThread_) {
    calcImage();
    ((SimulatorThread*)parent())->imageThreadDone();

    condition_.wait(&mutex_);
  }

  mutex_.unlock();
}

void ConcentImageThread::updateImage() {
  mutex_.lock();
  condition_.wakeOne();
  mutex_.unlock();
}

void ConcentImageThread::stopThread() {
  mutex_.lock();
  endThread_ = true;
  condition_.wakeOne();
  mutex_.unlock();
  wait();
}

void ConcentImageThread::calcImage() {
  const SimState &simState = simulator_.simulatedState();
  int nProds = simState.nProducts();

  for (int j = 0; j < size_.height(); ++j) {
    for (int i = 0; i < size_.width(); ++i) {
      int r=0, g=0, b=0;

      for (int k = 0; k < nProds; ++k) {
        double concentFact = simState.product(k)(i,j) / factor_;

        QRgb color = ColorsConfig::colorsMorp[k % ColorsConfig::nColors];

        r += qRed(color) * concentFact;
        g += qGreen(color) * concentFact;
        b += qBlue(color) * concentFact;
      }

      if (r>255)
        r = 255;
      if (g>255)
        g = 255;
      if (b>255)
        b = 255;

      ((QRgb*) image_.scanLine(j))[i] = qRgb(r, g, b);
    }
  }
}

}