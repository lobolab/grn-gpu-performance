// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#pragma once


#include <QThread>
#include <QMutex>
#include <QWaitCondition>
#include <QImage>

namespace LoboLab {

class SimulatorThread;
class Simulator;

class ConcentImageThread : public QThread {
  Q_OBJECT

 public:
  ConcentImageThread(Simulator &simulator, double concentImageFactor,
                     SimulatorThread *parent = 0);
  ~ConcentImageThread();

  inline double factor() const { return factor_; }
  inline void setFactor(double factor) { factor_ = factor; }

  void updateImage(); // Uses its own thread
  void calcImage(); // Uses caller thread

  void stopThread();

  inline const QImage &image() const { return image_; }

 protected:
  void run();

 private:

  QMutex mutex_;
  QWaitCondition condition_;

  Simulator &simulator_;
  QSize size_;
  QImage image_;
  double factor_;

  bool endThread_;
};

} // namespace LoboLab
