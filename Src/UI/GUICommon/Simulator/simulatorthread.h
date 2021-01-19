// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#pragma once

#include <QMutex>
#include <QSize>
#include <QThread>
#include <QWaitCondition>
#include <QImage>

#include "cellimagethread.h"
#include "concentimagethread.h"

namespace LoboLab {

class Simulator;

class SimulatorThread : public QThread {
  friend class CellImageThread;
  friend class ConcentImageThread;
  friend class GapImageThread;
  Q_OBJECT

 public:
  SimulatorThread(Simulator *simulator, double concentImageFactor, 
                  QObject *parent = 0);
  ~SimulatorThread();

  inline const QImage &getCellImage() const {
    return cellImageThread_->image();
  }

  inline const QImage &getConcentImage() const {
    return concentImageThread_->image();
  }

  inline double getConcentImageFactor() const {
    return concentImageThread_->factor();
  }

  inline void setConcentImageFactor(double factor) {
    concentImageThread_->setFactor(factor);
  }

  inline void setAutoUpdateImages(bool v) {
    autoUpdateImages_ = v;
  }

  inline int getNumStepsPerCycle() const {
    return nStepsPerCycle_;
  }

  inline void setNumStepsPerCycle(int v) {
    nStepsPerCycle_ = v;
  }

  void stopThread();
  void simNextCycle();
  void updateImages(); // blocking function: the calculation is done in the
  // caller's thread

 signals:
  void simCycleDone(double change);

 protected:
  void run();

 private:
  void imageThreadDone();

  QMutex mutex_;
  QWaitCondition thisThreadCond_;
  QWaitCondition imageThreadsCond_;

  CellImageThread *cellImageThread_;
  ConcentImageThread *concentImageThread_;

  Simulator *simulator_;

  bool autoUpdateImages_;
  bool endThread_;
  int nImageThreadsDone_;
  int nStepsPerCycle_;
};

} // namespace LoboLab
