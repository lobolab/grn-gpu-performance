// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#pragma once

#include <QGraphicsScene>

namespace LoboLab {

class ImageSaver {
 public:
  ~ImageSaver();

  static void saveScene(QGraphicsScene &scene, bool reversed, 
                        QWidget *parent = 0);
  static void saveScene(QGraphicsScene &scene, const QString &fileName, 
                        bool reversed);

  static void savePixmap(const QPixmap &pixmap, QWidget *parent = 0);
  static void savePixmap(const QPixmap &pixmap, const QString &fileName);

 private:
  static void crop(QImage &img, QRgb blankPixel);

  ImageSaver();
};

} // namespace LoboLab
