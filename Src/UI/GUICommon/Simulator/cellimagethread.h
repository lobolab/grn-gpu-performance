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

class CellImageThread : public QThread {
  Q_OBJECT

 public:
  CellImageThread(const Simulator &sim, SimulatorThread *parent);
  ~CellImageThread();

  void updateImage(); // Uses its own thread
  void calcImage(); // Uses caller thread

  void stopThread();

  inline const QImage &image() const { return image_; }

 protected:
  void run();

 private:

  QMutex mutex_;
  QWaitCondition condition_;

  const Simulator &simulator_;
  QSize size_;
  QImage image_;

  bool endThread_;
};

} // namespace LoboLab
