// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#include "simulatorthread.h"

#include "Simulator/simulator.h"
#include "Simulator/simstate.h"
#include <QDateTime>

namespace LoboLab {

SimulatorThread::SimulatorThread(Simulator *sim, double concentImageFactor,
                                 QObject *parent)
  : QThread(parent), 
    simulator_(sim), 
    autoUpdateImages_(true),
    nStepsPerCycle_(1) {
  cellImageThread_ = new CellImageThread(*simulator_, this);
  concentImageThread_ = new ConcentImageThread(*simulator_, concentImageFactor, 
                                               this);
}


SimulatorThread::~SimulatorThread() {
  stopThread();
  delete cellImageThread_;
  delete concentImageThread_;
}

void SimulatorThread::run() {
  mutex_.lock();
  endThread_ = false;
  nImageThreadsDone_ = 0;

  cellImageThread_->start();
  concentImageThread_->start();

  thisThreadCond_.wait(&mutex_); // waiting for image threads

  while (!endThread_) {
    double change = simulator_->simulate(nStepsPerCycle_);

    mutex_.unlock();
    cellImageThread_->updateImage();
    concentImageThread_->updateImage();
    mutex_.lock();
    thisThreadCond_.wait(&mutex_); // waiting for updateImage

    emit simCycleDone(change);

    if (!endThread_)
      thisThreadCond_.wait(&mutex_); //waiting for nextCycle
  }

  cellImageThread_->stopThread();
  concentImageThread_->stopThread();

  mutex_.unlock();
}

void SimulatorThread::stopThread() {
  mutex_.lock();
  endThread_ = true;
  thisThreadCond_.wakeOne();
  mutex_.unlock();
  wait();

  cellImageThread_->stopThread();
  concentImageThread_->stopThread();
}

void SimulatorThread::simNextCycle() {
  mutex_.lock();
  thisThreadCond_.wakeOne();
  mutex_.unlock();
}

void SimulatorThread::imageThreadDone() {
  mutex_.lock();

  if (nImageThreadsDone_ == 1) {
    nImageThreadsDone_ = 0;
    thisThreadCond_.wakeOne();
  } else
    ++nImageThreadsDone_;

  mutex_.unlock();
}

void SimulatorThread::updateImages() {
  cellImageThread_->calcImage();
  concentImageThread_->calcImage();
}

}