// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#include "simparams.h"
#include "DB/db.h"

namespace LoboLab {

SimParams::SimParams(int id, DB *db)
  : ed("SimParams", id, db) {
  load();
}

SimParams::~SimParams() {
}

SimParams::SimParams(const SimParams &source, bool maintainId)
  : ed(source.ed, maintainId) {
  copy(source);
}

SimParams &SimParams::operator=(const SimParams &source) {
  ed = source.ed;
  copy(source);

  return *this;
}

void SimParams::copy(const SimParams &source) {
  name = source.name;
  size = source.size;
  dt = source.dt;
  dx = source.dx;
  NumSimSteps = source.NumSimSteps;
  nSims = source.nSims;
  nInputMorphogens = source.nInputMorphogens;
  nTargetMorphogens = source.nTargetMorphogens;
  nMaxProducts = source.nMaxProducts;
  nMaxLinks = source.nMaxLinks;
  minConc = source.minConc;
  minConcChange = source.minConcChange;
  distErrorThreshold = source.distErrorThreshold;
  errorPrecision = source.errorPrecision;
  nCPUSlaves = source.nCPUSlaves;
  nGPUSlaves = source.nGPUSlaves;
  nGPUThreads = source.nGPUThreads;
  nTestIndividuals = source.nTestIndividuals;
  multiGPU = source.multiGPU;
}

// Persistence methods

void SimParams::load() {
  name = ed.loadValue(FName).toString();
  size.setWidth(ed.loadValue(FSizeX).toInt());
  size.setHeight(ed.loadValue(FSizeY).toInt());
  dt = ed.loadValue(FDt).toDouble();
  dx = ed.loadValue(FDx).toDouble();
  NumSimSteps = ed.loadValue(FNumSimSteps).toInt();
  nSims = ed.loadValue(FNumSims).toInt();
  nInputMorphogens = ed.loadValue(FNumInputMorphogens).toInt();
  nTargetMorphogens = ed.loadValue(FNumTargetMorphogens).toInt();
  nMaxProducts = ed.loadValue(FNumMaxProducts).toInt();
  nMaxLinks = ed.loadValue(FNumMaxLinks).toInt();
  minConc = ed.loadValue(FMinConc).toDouble();
  minConcChange = ed.loadValue(FMinConcChange).toDouble();
  distErrorThreshold = ed.loadValue(FDistErrorThreshold).toDouble();
  errorPrecision = ed.loadValue(FErrorPrecision).toInt();
  nCPUSlaves = ed.loadValue(FNumCPUSlaves).toInt();
  nGPUSlaves = ed.loadValue(FNumGPUSlaves).toInt();
  nGPUThreads = ed.loadValue(FNumGPUThreads).toInt();
  nTestIndividuals = ed.loadValue(FNumTestIndividuals).toInt();
  multiGPU = ed.loadValue(FMultiGPU).toBool();
  ed.loadFinished();
}

int SimParams::submit(DB *db) {
  QHash<QString, QVariant> values;
  values.insert("Name", name);
  values.insert("DimX", size.width());
  values.insert("DimY", size.height());
  values.insert("Dt", dt);
  values.insert("Dx", dx);
  values.insert("NumSimSteps", NumSimSteps);
  values.insert("NumSims", nSims);
  values.insert("NumInputMorphogens", nInputMorphogens);
  values.insert("NumTargetMorphogens", nTargetMorphogens);
  values.insert("NumMaxProducts", nMaxProducts);
  values.insert("NumMaxLinks", nMaxLinks);
  values.insert("MinConc", minConc);
  values.insert("MinConcChange", minConcChange);
  values.insert("DistErrorThreshold", distErrorThreshold);
  values.insert("ErrorPrecision", errorPrecision);
  values.insert("NumCPUSlaves", nCPUSlaves);
  values.insert("NumGPUSlaves", nGPUSlaves);
  values.insert("NumGPUThreads", nGPUThreads);
  values.insert("NumTestIndividuals", nTestIndividuals);
  values.insert("MultiGPU", multiGPU);

  return ed.submit(db, values);
}

bool SimParams::erase() {
  return ed.erase();
}
}
