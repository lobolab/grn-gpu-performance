// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#include "errorcalculatormultithread.h"
#include "Search/individual.h"
#include "Search/search.h"
#include "Search/searchparams.h"
#include "Common/log.h"
#include "Common/mathalgo.h"
#include "GPUSimulator/cevaluatormorphogens.h"
#include "Search/evaluatormorphogens.h"
#include "Simulator/simparams.h"

#include <iostream>
#include <time.h>

namespace LoboLab {

ErrorCalculatorMultiThread::ErrorCalculatorMultiThread(int nDemes, 
  int nCPUThreads, int nGPUThreads, const Search &search)
    : ErrorCalculator(), 
      nDemes_(nDemes),
      nIndQueuedDeme_(nDemes_, 0),
      nIndPendDeme_(nDemes_, 0),
      nToInitialize_(nCPUThreads + nGPUThreads) {


  mutex_.lock();

  for (int i = 0; i < nCPUThreads; ++i) {

    
    CalculatorThread *thread = new CalculatorThreadCPU(search, this);
    calculatorThreads_.append(thread);
    thread->start();
  }

  for (int i = 0; i < nGPUThreads; ++i) {

    CalculatorThread *thread = new CalculatorThreadGPU(i, search, this);
    calculatorThreads_.append(thread);
    thread->start();
  }

  // Wait for threads to finish initialization
  parentCondition_.wait(&mutex_);
  mutex_.unlock();
}

ErrorCalculatorMultiThread::~ErrorCalculatorMultiThread(void) {
  for (int i = 0; i < calculatorThreads_.size(); ++i) {
    CalculatorThread *thread = calculatorThreads_.at(i);
    thread->stopThread();
    delete thread;
  }
}


int ErrorCalculatorMultiThread::testPerformance(const Search &search, int numSimulations) {
  // 4-product model
  QString modelStr("(87889 0 0.614236 0.459672 0.912483|132569 0 0.842144 0.325605 0.630702|127105 0 0.217604 0.916419 0.895754|1 -2 0.828329 0.916419 0|111594 0 0.416697 0.325605 0.413702|64653 0 0.0174728 0.916419 0|108926 0 0.733034 0.826109 0|67959 0 0.812301 0.916419 0|0 -1 0.985517 0.00487246 0.516741|31044 0 0.937029 0.325605 0.610245|95513 0 0.986953 0.369568 0|107652 0 0.0174728 0.774791 0|136389 0 0.0174728 0.774791 0*64653 87889 0.356057 3.29451 1|132569 132569 0.605912 -4.91105 1|64653 127105 0.101022 -1.10124 1|0 1 0.0162271 3.14285 0|1 1 0.0162271 2.03104 0|107652 1 0.700604 -1.58752 1|111594 111594 0.23574 2.13902 0|0 64653 0.0723652 4.98268 0|108926 108926 0.537059 4.1708 1|64653 67959 0.510806 3.3672 1|64653 0 0.928964 2.75749 0|0 31044 0.6394 -2.10032 1|64653 31044 0.354384 2.20191 0|64653 95513 0.132222 1.8407 0|107652 95513 0.608596 -4.13056 1|1 107652 0.0694555 3.95009 0|64653 107652 0.94844 -4.28828 1|136389 111594 0.605417 -2.11831 1|87889 136389 0.320511 -1.59639 1)");
  double expectedError = 0.048701011019999997098129;
  
  int nMorphogens = search.simParams()->nInputMorphogens +
    search.simParams()->nTargetMorphogens;

  QList<Individual*> individuals;
  for (int i = 0; i < numSimulations; ++i) {
    Model *model = new Model();
    model->loadFromString(modelStr);
    Individual *ind = new Individual(model, nMorphogens);
    individuals.append(ind);
  }

  process(0, individuals);
  waitForAnyDeme();

  int ret = 0;
  for (int i = 0; i < numSimulations; ++i) {
    if (double errorDif = qAbs(individuals[i]->error() - expectedError)) {
      Log::write() << "WRONG ERROR: ErrorDif = " << errorDif << ", Obtained Error = " << individuals[i]->error() << ", Expected Error = " << expectedError << endl;
      ret = -11;
    }

    delete individuals[i];
  }
  Log::write() << "ErrorCalculatorMultiThread::testPerformance: done. All errors checked." << endl;

  return ret;
}



void ErrorCalculatorMultiThread::process(int iDeme,
                                        const QList<Individual*> &individuals) {
  mutex_.lock();

  pendDemeQueue_.enqueue(iDeme);
  pendIndQueue_.append(individuals);
  nIndQueuedDeme_[iDeme] = individuals.size();
  nIndPendDeme_[iDeme] = individuals.size();

  if (!waitingThreads_.isEmpty())
    waitingThreads_.dequeue()->wake();
  
  mutex_.unlock();
}

int ErrorCalculatorMultiThread::waitForAnyDeme() {
  mutex_.lock();

  while (readyDemeQueue_.isEmpty())
    parentCondition_.wait(&mutex_);

  int iDemeReady = readyDemeQueue_.dequeue();

  mutex_.unlock();

  return iDemeReady;
}

// class CalculatorThread

ErrorCalculatorMultiThread::CalculatorThread::CalculatorThread(
    const Search &search,
    ErrorCalculatorMultiThread *p)
  : parent_(p) {
  evalmorph_ = new EvaluatorMorphogens(search);
}

ErrorCalculatorMultiThread::CalculatorThread::~CalculatorThread() {
  stopThread();
  delete evalmorph_;
}

// Waits for all threads to initialize, so memory allocations are finished
// before launching kernels.
// Children can reimplemnt this virtual function, but they need to call this
// parent function when finish initializing.
void ErrorCalculatorMultiThread::CalculatorThread::initialize() {
  parent_->mutex_.lock();

  if (--parent_->nToInitialize_ > 0) {
    parent_->parentCondition_.wait(&parent_->mutex_);
  }
  else {
    parent_->parentCondition_.wakeAll();
  }

  parent_->mutex_.unlock();
}

void ErrorCalculatorMultiThread::CalculatorThread::run() {
  initialize();
  setPriority(LowestPriority);

  parent_->mutex_.lock();
  endThread_ = false;

  while (!endThread_) {
    if (parent_->pendIndQueue_.isEmpty())
      wait();
    else {
      processNextIndividual();
      if (!parent_->waitingThreads_.isEmpty()) {
        parent_->waitingThreads_.dequeue()->wake();
        wait();
      }
    }
  }

  parent_->mutex_.unlock();
}

void ErrorCalculatorMultiThread::CalculatorThread::processNextIndividual() {
  Individual *nextInd = parent_->pendIndQueue_.dequeue();
  
  int iDeme = parent_->pendDemeQueue_.head();

  if (--parent_->nIndQueuedDeme_[iDeme] == 0) // last individual in queue from deme
    parent_->pendDemeQueue_.dequeue();

  if (!parent_->pendIndQueue_.isEmpty() && !parent_->waitingThreads_.isEmpty())
    parent_->waitingThreads_.dequeue()->wake(); // more individuals to process

  parent_->mutex_.unlock();

  double error, simTime;
  calcError(*nextInd->model(), nextInd->parentError(), &error, &simTime);
  nextInd->setError(error);
  nextInd->setSimTime(simTime);
  
  parent_->mutex_.lock();

  if (--parent_->nIndPendDeme_[iDeme] == 0) { // last individual processed in deme
    parent_->readyDemeQueue_.enqueue(iDeme);
    parent_->parentCondition_.wakeOne();
  }
}

void ErrorCalculatorMultiThread::CalculatorThread::stopThread() {
  parent_->mutex_.lock();
  endThread_ = true;
  wake();
  parent_->mutex_.unlock();
  QThread::wait();
}

// class CalculatorThreadCPU

//Same constructor as CalculatorThread
ErrorCalculatorMultiThread::CalculatorThreadCPU::
CalculatorThreadCPU(const Search &search, ErrorCalculatorMultiThread *p) 
  : CalculatorThread(search, p) {
  
}

ErrorCalculatorMultiThread::CalculatorThreadCPU::~CalculatorThreadCPU(void) {
  
}

//Copied from ErrorCalculatorMultiThread CPU Version
void ErrorCalculatorMultiThread::CalculatorThreadCPU::
calcError(const Model &model, double maxError, double *error, double *simTime) {

  timer_.start();
  *error = evalmorph_->evaluate(model, maxError);
  *simTime = timer_.elapsed() / 1000.0;
}

// class CalculatorThreadGPU

ErrorCalculatorMultiThread::CalculatorThreadGPU::
CalculatorThreadGPU(int threadId, const Search &search, ErrorCalculatorMultiThread *p)
  : CalculatorThread(search, p), threadId_(threadId), search_(search), cevalmorph_(NULL) {
}

ErrorCalculatorMultiThread::CalculatorThreadGPU::~CalculatorThreadGPU(void) {
  delete cevalmorph_;
}

void ErrorCalculatorMultiThread::CalculatorThreadGPU::initialize() {
  cevalmorph_ = new CEvaluatorMorphogens(threadId_, search_);
  CalculatorThread::initialize();
}

//Copied from ErrorCalculatorMultiThread GPU Version
void ErrorCalculatorMultiThread::CalculatorThreadGPU::
calcError(const Model &model, double maxError, double *error, double *simTime) {

  timer_.start();
  *error = cevalmorph_->evaluate(model, maxError);
  *simTime = timer_.elapsed() / 1000.0;
}

}
