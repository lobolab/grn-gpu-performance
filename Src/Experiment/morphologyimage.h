// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#pragma once

#include "DB/dbelementdata.h"
#include <QByteArray>
#include <QPixmap>

namespace LoboLab {

class Experiment; 

class MorphologyImage : public DBElement {

  public:
    MorphologyImage(int id, DB *db);
    virtual ~MorphologyImage();
    QByteArray format();
    const QImage &image() const;
    const QPixmap &pixmap() const;

    inline const QString &name() const { return name_; }
    bool operator==(const MorphologyImage& other) const;
    inline bool operator!=(const MorphologyImage& other) const {
      return !(*this == other);
    }

  protected:
    inline virtual int id() const { return ed_.id(); };
    virtual int submit(DB *db);
    virtual bool erase();

  private:

    MorphologyImage & operator=(const MorphologyImage &source);
    
    void load();

    QByteArray data_;
    QString name_;

    mutable QImage *image_;
    mutable QPixmap *pixmap_;

    DBElementData ed_;

    // Persistence fields
  public:
    enum {
      FName = 1,
      FData = 2
    };
};

} // namespace LoboLab
