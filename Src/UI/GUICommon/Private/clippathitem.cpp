// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#include "clippathitem.h"

#include <QPen>

namespace LoboLab {

ClipPathItem::ClipPathItem(const QPainterPath &c, QGraphicsItem *child,
                           QGraphicsItem *parent)
  : QGraphicsPathItem(parent), clip_(c), bRect_(c.boundingRect()) {
  child->setParentItem(this);
  setFlag(ItemClipsChildrenToShape);
  setPen(Qt::NoPen);
  setBrush(Qt::NoBrush);
}

ClipPathItem::~ClipPathItem() {
}

QRectF ClipPathItem::boundingRect() const {
  return bRect_;
}

QPainterPath ClipPathItem::shape() const {
  return clip_;
}


}