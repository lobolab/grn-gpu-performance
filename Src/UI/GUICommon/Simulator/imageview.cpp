// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#include "imageview.h"
#include "UI/GUICommon/imagesaver.h"

#include <QPainter>
#include <QStyle>
#include <QMouseEvent>
#include <QAction>

namespace LoboLab {

ImageView::ImageView(bool showSelecLin, QWidget *parent)
    : QWidget(parent), 
      useSelecLines_(showSelecLin), 
      lineHpos_(0.51), lineVpos_(0.51),
      movingLineH_(false), movingLineV_(false), 
      drawingRect_(false), 
      rectColor_(Qt::white),
      drawBorder_(true) {
  setStyleSheet("QMenu {background-color:#F0F0F0;} * {background-color:transparent;}");
  QSizePolicy sp(QSizePolicy::Preferred, QSizePolicy::Preferred);
  sp.setHeightForWidth(true);
  setSizePolicy(sp);

  tempImage_ = new QImage(1, 1, QImage::Format_RGB32);
  image_ = tempImage_;

  QAction *saveImageAction_ = new QAction(tr("Save image..."), this);
  saveImageAction_->setIcon(QIcon(":/Images/famfamfam_silk_icons/picture_save.png"));
  connect(saveImageAction_, SIGNAL(triggered()), this, SLOT(saveImage()));

  setContextMenuPolicy(Qt::ActionsContextMenu);
  addAction(saveImageAction_); // add the action to the context menu
}

ImageView::~ImageView() {
  delete tempImage_;
}
void ImageView::setRectColor(const QColor &color) {
  rectColor_ = color;
}

void ImageView::setImage(const QImage *img) {
  image_ = img;
  setMinimumSize(image_->size());
  update();
}

void ImageView::paintEvent(QPaintEvent *event) {
  if (image_) {
    QStyle *style = QWidget::style();

    QPainter painter(this);
    QRect cr = contentsRect();
    float minS = MathAlgo::min(cr.width(), cr.height());

    int exactS = image_->width() * floor(minS / image_->width());
    
    QRect cr2((cr.width()-exactS)/2, (cr.width()-exactS)/2, exactS, exactS);
        
    style->drawItemPixmap(&painter, cr2, Qt::AlignCenter | Qt::AlignCenter,
                          QPixmap::fromImage(image_->scaled(cr2.size(), 
                                                            Qt::IgnoreAspectRatio,
                                                            Qt::FastTransformation),
                                             Qt::NoOpaqueDetection));
    
    if (useSelecLines_) {
      painter.setPen(Qt::lightGray);
      painter.drawLine(cr2.x(), cr2.y() + cr2.height()*lineHpos_,
                       cr2.x() + cr2.width() - 1, cr2.y() + cr2.height()*lineHpos_);

      painter.setPen(Qt::gray);
      painter.drawLine(cr2.x() + cr2.width()*lineVpos_, cr2.y(),
                       cr2.x() + cr2.width()*lineVpos_, cr2.y() + cr2.height() - 1);
    }

    if (drawingRect_) {
      painter.setPen(rectColor_);
      painter.drawRect(QRect(lastPoint_, currentPoint_));
    }
    
    if (drawBorder_) {
      painter.setPen(Qt::lightGray);
      painter.drawRect(cr.adjusted(0,0,-1,-1));
    }
  }
}

void ImageView::mousePressEvent(QMouseEvent *event) {
  lastPoint_ = currentPoint_ = event->pos();

  if (useSelecLines_ && event->button() == Qt::LeftButton) {
    int margin = 10;
    QRect cr = contentsRect();

    if (currentPoint_.x() < lineVpos_ * cr.width() + margin &&
        currentPoint_.x() > lineVpos_ * cr.width() - margin)
      movingLineV_ = true;

    if (currentPoint_.y() < lineHpos_ * cr.height() + margin &&
        currentPoint_.y() > lineHpos_ * cr.height() - margin)
      movingLineH_ = true;

    if (!movingLineV_ && !movingLineH_)
      drawingRect_ = true;
  } else
    drawingRect_ = true;

  event->accept();
}

void ImageView::mouseMoveEvent(QMouseEvent *event) {
  currentPoint_ = event->pos();

  if (movingLineH_) {
    QRect cr = contentsRect();
    lineHpos_ = (float)currentPoint_.y()/cr.height();

    if (lineHpos_ < 0)
      lineHpos_ = 0;
    else if (lineHpos_ > 1)
      lineHpos_ = 1;

    update();
    emit lineHmoved();
  }

  if (movingLineV_) {
    QRect cr = contentsRect();
    lineVpos_ = (float)currentPoint_.x()/cr.width();

    if (lineVpos_ < 0)
      lineVpos_ = 0;
    else if (lineVpos_ > 1)
      lineVpos_ = 1;

    update();
    emit lineVmoved();
  }

  if (drawingRect_) {
    update();
  }

  event->accept();
}

void ImageView::wheelEvent(QWheelEvent *event) {
  emit mouseWheel(event);
  event->accept();
}

void ImageView::mouseReleaseEvent(QMouseEvent *event) {
  if (movingLineH_)
    movingLineH_ = false;
  if (movingLineV_)
    movingLineV_ = false;

  if (drawingRect_) {
    QRect rect = QRect(lastPoint_, event->pos());
    emit rectDrawed(QRect(rect.topLeft(), rect.bottomRight()),
                    event->button());
    drawingRect_ = false;
    update();
  }
}

// private slot
void ImageView::saveImage() {
  ImageSaver::savePixmap(createPixmap(), parentWidget());
}

QPixmap ImageView::createPixmap() {
  QPixmap pixmap(size());
  pixmap.fill(Qt::transparent);

  bool prevUseSelecLines = useSelecLines_;
  useSelecLines_ = false;
  drawBorder_ = false;
  render(&pixmap);
  drawBorder_ = true;
  useSelecLines_ = prevUseSelecLines;

  return pixmap;
}

}