// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#include "experiment.h"
#include "morphologyimage.h"
#include <string>
#include <iostream>


namespace LoboLab {

Experiment::Experiment()
  : ed_("Experiment") {
}

Experiment::Experiment(int id, DB *db) // loadImages is added
  : ed_("Experiment", id, db) {
  load();
}

Experiment::~Experiment() {
  deleteAll();
}

Experiment::Experiment(const Experiment &source, bool maintainId)
  : ed_(source.ed_, maintainId) {
  copy(source, maintainId);
}

Experiment &Experiment::operator=(const Experiment &source) {
  deleteAll();
  ed_ = source.ed_;
  copy(source, true);

  return *this;
}

Experiment &Experiment::copy(const Experiment &source) {
  removeAndDeleteAll();
  copy(source, false);

  return *this;
}

void Experiment::copy(const Experiment &source, bool maintainId) {
  name_ = source.name_;
  inputMorphology_ = new MorphologyImage(*source.inputMorphology_);
  outputMorphology_ = new MorphologyImage(*source.outputMorphology_);
}

void Experiment::deleteAll() {

  delete inputMorphology_;
  inputMorphology_ = NULL;
  delete outputMorphology_;
  outputMorphology_ = NULL;

}

void Experiment::setInputMorphology(MorphologyImage *i) {
  delete inputMorphology_;
  inputMorphology_ = i;
}

void Experiment::setOutputMorphology(MorphologyImage *o) {
  delete outputMorphology_;
  outputMorphology_ = o;
}

void Experiment::removeAndDeleteAll() {

  delete inputMorphology_;
  inputMorphology_ = NULL;
  delete outputMorphology_;
  outputMorphology_ = NULL;

}

bool Experiment::operator==(const Experiment& other) const {
  bool equal = name_ == other.name_ &&
    inputMorphology_->name() == other.inputMorphology_->name() &&
    outputMorphology_->name() == other.outputMorphology_->name();
  return equal;
}

// Persistence methods

void Experiment::load() {
  name_ = ed_.loadValue(FName).toString();
  inputMorphology_ = new MorphologyImage(ed_.loadValue(FInputMorphology).toInt(), ed_.db());
  outputMorphology_ = new MorphologyImage(ed_.loadValue(FOutputMorphology).toInt(), ed_.db());

  ed_.loadFinished();
}

int Experiment::submit(DB *db) {
  QPair<QString, DBElement*> refMember;

  QHash<QString, QVariant> values;
  values.insert("Name", name_);

  QHash<QString, DBElement*> members;
  members.insert("InputMorphology", inputMorphology_);
  members.insert("OutputMorphology", outputMorphology_);

  return ed_.submit(db, refMember, members, values);
}

bool Experiment::erase() {
  QList<DBElement*> members;

  return ed_.erase(members);
}

}