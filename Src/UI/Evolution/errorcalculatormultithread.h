// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#pragma once

#include "GPUSimulator/cevaluatormorphogens.h"
#include "Search/errorcalculator.h"
#include "Search/search.h"

#include <QList>
#include <QThread>
#include <QMutex>
#include <QWaitCondition>
#include <QQueue>
#include <QVector>
#include <QElapsedTimer>

namespace LoboLab {

class EvaluatorMorphogens;
class Model;
class DB;
class SimParams;
class CEvaluatorMorphogens;

class ErrorCalculatorMultiThread : public ErrorCalculator {

 public:
  ErrorCalculatorMultiThread(int nDemes, int nCPUThreads, int nGPUThreads, 
                              const Search &search);
  virtual ~ErrorCalculatorMultiThread(void);
 
  void process(int iDeme, const QList<Individual*> &individuals);
  int waitForAnyDeme();

  int testPerformance(const Search &search, int numSimulations);


 private:

   class CalculatorThread : public QThread {
   public:
     CalculatorThread(const Search &search, ErrorCalculatorMultiThread *parent);
     ~CalculatorThread(void);

     virtual void initialize();
     void stopThread();
     inline void wake() { waitCondition_.wakeOne();  };
     inline void wait() { parent_->waitingThreads_.enqueue(this);  waitCondition_.wait(&parent_->mutex_); };;

   protected:
     void run();
     EvaluatorMorphogens *evalmorph_;
     
     QElapsedTimer timer_;

   private:

    void processNextIndividual();
    virtual void calcError(const Model &model, double maxError, double *error,
                           double *simTime) = 0;

    ErrorCalculatorMultiThread *parent_;
    Individual* individual_;

    bool endThread_;
    QWaitCondition waitCondition_;
  };
   
   class CalculatorThreadCPU : public CalculatorThread {
   public:
     CalculatorThreadCPU(const Search &search,
                         ErrorCalculatorMultiThread *parent);
     ~CalculatorThreadCPU(void);

   private:
     void calcError(const Model &model, double maxError, double *error,
                    double *simTime);

   };

   class CalculatorThreadGPU : public CalculatorThread {
   public:
     CalculatorThreadGPU(int threadId, const Search &search,
                         ErrorCalculatorMultiThread *parent);
     ~CalculatorThreadGPU(void);

     virtual void initialize();

   private:
     void calcError(const Model &model, double maxError, double *error,
                    double *simTime);
     
     int threadId_;
     const Search &search_;
     CEvaluatorMorphogens *cevalmorph_;
   };

  int nDemes_;
  QQueue<int> pendDemeQueue_;
  QQueue<int> readyDemeQueue_;
  QQueue<Individual*> pendIndQueue_;
  QVector<int> nIndQueuedDeme_;
  QVector<int> nIndPendDeme_;
  
  QMutex mutex_;
  QWaitCondition parentCondition_;
  int nToInitialize_;

  QList<CalculatorThread*> calculatorThreads_;
  QQueue<CalculatorThread*> waitingThreads_;
};


} // namespace LoboLab
