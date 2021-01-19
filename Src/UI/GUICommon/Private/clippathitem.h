// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#pragma once

#include <QGraphicsPathItem>

namespace LoboLab {

class ClipPathItem : public QGraphicsPathItem {
 public:
  ClipPathItem(const QPainterPath &clip, QGraphicsItem *child,
               QGraphicsItem *parent = NULL);
  virtual ~ClipPathItem();

  virtual QRectF boundingRect() const;
  virtual QPainterPath shape() const;

 private:
  QPainterPath clip_;
  QRectF bRect_;
};

} // namespace LoboLab
