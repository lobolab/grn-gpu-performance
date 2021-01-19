// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#pragma once

#include "Simulator/simstate.h"

#include <QWidget>
#include <QColor>
#include <QImage>
#include <QPoint>

namespace LoboLab {

class ImageView : public QWidget {
  Q_OBJECT

 public:
  ImageView(bool useSelecLines = false, QWidget *parent = NULL);
  ~ImageView();

  void setRectColor(const QColor &color);
  void setImage(const QImage *image);

  inline float lineHpos() const { return lineHpos_; }
  inline float lineVpos() const { return lineVpos_; }

  inline void setLineHpos(float pos) { lineHpos_ = pos; }
  inline void setLineVpos(float pos) { lineVpos_ = pos; }

  // Keep aspect ratio
  virtual int heightForWidth(int w) const { return w; }

  QPixmap createPixmap();

 signals:
  void rectDrawed(const QRect &rect, Qt::MouseButton button);
  void mouseWheel(QWheelEvent *event);
  void lineHmoved();
  void lineVmoved();

 private slots:
  void saveImage();

 protected:
  virtual void paintEvent(QPaintEvent *event);
  virtual void mousePressEvent(QMouseEvent *event);
  virtual void mouseMoveEvent(QMouseEvent *event);
  virtual void mouseReleaseEvent(QMouseEvent *event);
  virtual void wheelEvent(QWheelEvent *event);

 private:
  bool useSelecLines_;
  double lineHpos_;
  double lineVpos_;
  bool movingLineH_;
  bool movingLineV_;
  bool drawingRect_;
  QPoint lastPoint_;
  QPoint currentPoint_;
  QImage *tempImage_;
  const QImage *image_;
  QColor rectColor_;
  bool drawBorder_;
};

} // namespace LoboLab
