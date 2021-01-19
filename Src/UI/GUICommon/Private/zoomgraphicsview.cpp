// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#include "zoomgraphicsview.h"
#include "UI/GUICommon/imagesaver.h"

#include <QWheelEvent>
#include <QAction>

namespace LoboLab {

ZoomGraphicsView::ZoomGraphicsView(QWidget *parent, bool keepFit, bool rev)
  : QGraphicsView(parent), 
    keepFitted_(keepFit), 
    reversed_(rev), 
    zoomMode_(false),
    fittingView_(false) {
  scene_ = new QGraphicsScene(this);
  scene_->setSceneRect(1, 1, 1, 1);
  setScene(scene_);
  setRenderHints(QPainter::Antialiasing | QPainter::TextAntialiasing |
                 QPainter::SmoothPixmapTransform);
  if (reversed_)
    scale(1, -1); // adjust the axis to grow from left to right (x)
  // and from down to top (y)

  //setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
  //setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);

  saveImageAction_ = new QAction(tr("Save image..."), this);
  saveImageAction_->setIcon(QIcon(":/Images/famfamfam_silk_icons/picture_save.png"));
  connect(saveImageAction_, SIGNAL(triggered()), this, SLOT(saveImage()));

  setContextMenuPolicy(Qt::ActionsContextMenu);
  addAction(saveImageAction_); // add the action to the context menu
}

ZoomGraphicsView::~ZoomGraphicsView() {
  delete scene_;
}

void ZoomGraphicsView::wheelEvent(QWheelEvent* event) {
  //Get the position of the mouse before scaling, in scene coords
  //QPointF pointBeforeScale(mapToScene(event->pos()));

  //Get the original screen centerpoint
  //QPointF screenCenter = GetCenter(); //CurrentCenterPoint; //(visRect.center());

  //Scale the view ie. do the zoom
  double scaleFactor = 1.15; //How fast we zoom
  if (zoomMode_) {
    if (event->delta() > 0) {
      //Zoom in
      scale(scaleFactor, scaleFactor);
    } else {
      //Zooming out
      scale(1.0 / scaleFactor, 1.0 / scaleFactor);
    }
  } else
    QGraphicsView::wheelEvent(event);

  //Get the position after scaling, in scene coords
  //QPointF pointAfterScale(mapToScene(event->pos()));

  //Get the offset of how the screen moved
  //QPointF offset = pointBeforeScale - pointAfterScale;

  //Adjust to the new center for correct zooming
  //QPointF newCenter = screenCenter + offset;
  //SetCenter(newCenter);
}

void ZoomGraphicsView::keyPressEvent(QKeyEvent *e) {
  if (e->key() == Qt::Key_Control)
    zoomMode_ = true;
}

void ZoomGraphicsView::keyReleaseEvent(QKeyEvent *e) {
  if (e->key() == Qt::Key_Control)
    zoomMode_ = false;
}

void ZoomGraphicsView::focusOutEvent(QFocusEvent *) {
  zoomMode_ = false;
}

// private slot
void ZoomGraphicsView::resizeEvent(QResizeEvent *event) {
  QGraphicsView::resizeEvent(event);
  if (!fittingView_ && !zoomMode_)
    fitView();
}

void ZoomGraphicsView::fitView() {
  if (!fittingView_) {
    fittingView_ = true;

    setTransform(QTransform());
    if (reversed_)
      scale(1, -1);

    const QRectF sceneRect = scene_->itemsBoundingRect(); // scene->sceneRect();
    scene_->setSceneRect(sceneRect);

    if (keepFitted_) {
      const QRectF viewRect = rect();
      if (sceneRect.width() > viewRect.width() ||
          sceneRect.height() > viewRect.height()) {
        fitInView(sceneRect, Qt::KeepAspectRatio);
      }
    }

    fittingView_ = false;
  }
}

// private slot
void ZoomGraphicsView::saveImage() {
  ImageSaver::saveScene(*scene_, reversed_, parentWidget());
}

}
