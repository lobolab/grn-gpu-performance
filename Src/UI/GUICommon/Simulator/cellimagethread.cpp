// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#include "cellimagethread.h"

#include "simulatorthread.h"

#include "Simulator/simulator.h"
#include "Simulator/simstate.h"

#include "../Private/colorsconfig.h"
#include <QImage>
#include "Common/log.h"

namespace LoboLab {

CellImageThread::CellImageThread(const Simulator &sim,
                                 SimulatorThread *parent)
  : QThread(parent), simulator_(sim), size_(simulator_.domainSize()),
    image_(size_, QImage::Format_ARGB32_Premultiplied) {
  calcImage();
}

CellImageThread::~CellImageThread() {
  stopThread();
}

void CellImageThread::run() {
  mutex_.lock();
  endThread_ = false;

  while (!endThread_) {
    calcImage();
    ((SimulatorThread*)parent())->imageThreadDone();

    condition_.wait(&mutex_);
  }

  mutex_.unlock();
}

void CellImageThread::updateImage() {
  mutex_.lock();
  condition_.wakeOne();
  mutex_.unlock();
}

void CellImageThread::stopThread() {
  mutex_.lock();
  endThread_ = true;
  condition_.wakeOne();
  mutex_.unlock();
  wait();
}

void CellImageThread::calcImage() {
  const SimState &simState = simulator_.simulatedState();
  int nCellTypeProducts = 3;

  for (int j = 0; j < size_.height(); ++j) {
    for (int i = 0; i < size_.width(); ++i) {
      QRgb color;
	    color = 0xFF000000;
	  
	  
		  int red=0, green=0, blue=0;
		  
		  switch (simState.nInputMorphogens() + simState.nTargetMorphogens()) {
	      case 3: {
		      blue = 255 * MathAlgo::min(1.0, simState.product(2)(i, j));
	      }
	      case 2: {
		      green = 255 * MathAlgo::min(1.0, simState.product(1)(i, j));
	      }
	      case 1: {
		      red = 255 * MathAlgo::min(1.0, simState.product(0)(i, j));
	      }
      }

		  color = qRgb(red, green, blue);

      ((QRgb*) image_.scanLine(j))[i] = color;
    }
  }
}

}