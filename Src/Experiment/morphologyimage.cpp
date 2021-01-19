// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#include "morphologyimage.h"
#include "experiment.h"
#include <QBuffer>
#include <QImageReader>
#include <QFile>

namespace LoboLab {

MorphologyImage::MorphologyImage(int id, DB *db)
  : ed_("MorphologyImage", id, db), image_(NULL), pixmap_(NULL) {
  load();
}

MorphologyImage::~MorphologyImage() {
  delete image_;
  delete pixmap_;
}

QByteArray MorphologyImage::format() {
  QBuffer buffer(&data_);
  buffer.open(QIODevice::ReadOnly);
  QImageReader imageReader(&buffer);

  return imageReader.format();
}

const QImage &MorphologyImage::image() const {
  if (!image_) {
    image_ = new QImage();
    bool ok = image_->loadFromData(data_);

    if (!ok) {
      *image_ = QImage(100, 100, QImage::Format_RGB32);
      image_->fill(Qt::red);
    }
  }

  return *image_;
}

const QPixmap &MorphologyImage::pixmap() const {
  if (!pixmap_) {
    pixmap_ = new QPixmap();
    bool ok = pixmap_->loadFromData(data_);

    if (!ok) {
      *pixmap_ = QPixmap(100, 100);
      pixmap_->fill(Qt::red);
    }
  }

  return *pixmap_;
}

bool MorphologyImage::operator==(const MorphologyImage& other) const {
  return data_ == other.data_ &&
    name_ == other.name_;
}

// Persistence methods
void MorphologyImage::load() {
  name_ = ed_.loadValue(FName).toString();
  data_ = ed_.loadValue(FData).toByteArray();

  ed_.loadFinished();
}

int MorphologyImage::submit(DB *db) {

  QHash<QString, QVariant> values;
  values.insert("Data", data_);
  values.insert("Name", name_);

  return ed_.submit(db, values);
}

bool MorphologyImage::erase() {
  return ed_.erase();
}

}