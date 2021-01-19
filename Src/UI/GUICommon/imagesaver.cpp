// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#include "imagesaver.h"
#include "Common/mathalgo.h"

#include <QGraphicsScene>
#include <QFileDialog>
#include <QPixmap>
#include <QWidget>
//#include <QSvgGenerator>
#include <QPainter>
#include <QMessageBox>
#include <QSettings>

namespace LoboLab {

ImageSaver::ImageSaver() {
}

ImageSaver::~ImageSaver() {
}

void ImageSaver::saveScene(QGraphicsScene &scene, bool reversed, 
                           QWidget *parent) {
  QSettings settings;
  settings.beginGroup("ImageSaver");
  QString lastDir = settings.value("lastDir").toString();

  QString fileName = QFileDialog::getSaveFileName(parent, "Save image",
                     lastDir, "Images (*.png *.jpg *.bmp)");

  if (!fileName.isEmpty()) {
    QString newDir = QFileInfo(fileName).absolutePath();
    if (newDir != lastDir)
      settings.setValue("lastDir", newDir);

    saveScene(scene, fileName, reversed);
  }
}

void ImageSaver::saveScene(QGraphicsScene &scene, const QString &fileName,
                           bool reversed) {
  QRectF rect = scene.itemsBoundingRect();
  const int maxSize = 25000;
  int maxSizeRect = MathAlgo::max(rect.width(), rect.height());
  double fac = ((double)maxSize)/maxSizeRect;
  if (fac > 2)
    fac = 2;
  
	int sizeX = MathAlgo::max(1.0, fac*rect.width());
	int sizeY = MathAlgo::max(1.0, fac*rect.height());

  QImage img(sizeX, sizeY, QImage::Format_ARGB32_Premultiplied);
  QRgb blankPixel;
  if (fileName.endsWith(".png")) { // pngs are saved with transparent background
    img.fill(Qt::transparent);
    blankPixel = 0;
  } else {
    img.fill(Qt::white);
    blankPixel =  0xFFFFFFFF;
  }

  QPainter painter(&img);
  painter.setRenderHints(QPainter::Antialiasing | QPainter::TextAntialiasing |
                         QPainter::SmoothPixmapTransform);
  scene.render(&painter, QRectF(), rect);
  painter.end();

  crop(img, blankPixel);

  bool ok;
  if(reversed)
    ok = img.mirrored().save(fileName);
  else
    ok = img.save(fileName);

  Q_ASSERT(ok);
  (void)ok;
}

void ImageSaver::crop(QImage &img, QRgb blankPixel) {
  QRect cropArea = img.rect();

  bool blankLine;

  // top
  blankLine = true;
  while (blankLine && cropArea.height() > 1) {
    for (int i = cropArea.left(); i <= cropArea.right() && blankLine; ++i)
      blankLine = img.pixel(i, cropArea.top()) == blankPixel;

    if (blankLine)
      cropArea.setTop(cropArea.top() + 1);
  }

  // bottom
  blankLine = true;
  while (blankLine && cropArea.height() > 1) {
    for (int i = cropArea.left(); i <= cropArea.right() && blankLine; ++i)
      blankLine = img.pixel(i, cropArea.bottom()) == blankPixel;

    if (blankLine)
      cropArea.setBottom(cropArea.bottom() - 1);
  }

  // left
  blankLine = true;
  while (blankLine && cropArea.width() > 1) {
    for (int i = cropArea.top(); i <= cropArea.bottom() && blankLine; ++i)
      blankLine = img.pixel(cropArea.left(), i) == blankPixel;

    if (blankLine)
      cropArea.setLeft(cropArea.left() + 1);
  }

  // right
  blankLine = true;
  while (blankLine && cropArea.width() > 1) {
    for (int i = cropArea.top(); i <= cropArea.bottom() && blankLine; ++i)
      blankLine = img.pixel(cropArea.right(), i) == blankPixel;

    if (blankLine)
      cropArea.setRight(cropArea.right() - 1);
  }

  if (cropArea != img.rect())
    img = img.copy(cropArea);
}

void ImageSaver::savePixmap(const QPixmap &pixmap, QWidget *parent) {
  QSettings settings;
  settings.beginGroup("ImageSaver");
  QString lastDir = settings.value("lastDir").toString();

  QString fileName = QFileDialog::getSaveFileName(parent, "Save image",
                     lastDir, "Images (*.png *.jpg *.bmp)");

  if (!fileName.isEmpty()) {
    QString newDir = QFileInfo(fileName).absolutePath();
    if (newDir != lastDir)
      settings.setValue("lastDir", newDir);

    savePixmap(pixmap, fileName);
  }
}

void ImageSaver::savePixmap(const QPixmap &pixmap, const QString &fileName) {
  QImage img = pixmap.toImage();
  crop(img, img.pixel(0,0));

  bool ok = img.save(fileName);
  Q_ASSERT(ok);
  (void)ok;
}

}
