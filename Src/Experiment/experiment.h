// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#pragma once

#include "DB/dbelementdata.h"

namespace LoboLab {

class MorphologyImage;

class Experiment : public DBElement {
  public:
    Experiment();
    Experiment(int id, DB *db);
    virtual ~Experiment();

    Experiment(const Experiment &source, bool maintainId = true);
    Experiment &operator=(const Experiment &source);
    Experiment &copy(const Experiment &source);

    inline const QString &name() const { return name_; }
    inline MorphologyImage *inputMorphology() const { return inputMorphology_; }
    inline MorphologyImage *outputMorphology() const { return outputMorphology_; }

    inline void setName(const QString &name) { name_ = name; }
    void setInputMorphology(MorphologyImage *inputMorphology);
    void setOutputMorphology(MorphologyImage *outputMorphology);


    bool operator==(const Experiment& other) const;
    inline bool operator!=(const Experiment& other) const {
      return !(*this == other);
    }

    inline virtual int id() const { return ed_.id(); }
    virtual int submit(DB *db);
    virtual bool erase();

  private:
    void copy(const Experiment &source, bool maintainId);
    void deleteAll();
    void removeAndDeleteAll();

    void load();

    QString name_;

    MorphologyImage *inputMorphology_;
    MorphologyImage *outputMorphology_;

    mutable DBElementData ed_;

  public:
    enum {
      FName = 1,
      FInputMorphology = 2,
      FOutputMorphology = 3
    };
};

} // namespace LoboLab
