// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#include "morphologyimagewidget.h"
#include "Experiment/morphologyimage.h"
#include "UI/GUICommon/imagesaver.h"

#include <QPen>
#include <QAction>

namespace LoboLab {

MorphologyImageWidget::MorphologyImageWidget(MorphologyImage *morphologyImage,
  QWidget *parent)
  : QLabel(parent), morphologyImage_(morphologyImage), pixmap_(morphologyImage->pixmap()) {
  setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Ignored);

  QAction *saveImageAction = new QAction(tr("Save image..."), this);
  saveImageAction->setIcon(QIcon(":/Images/famfamfam_silk_icons/picture_save.png"));
  connect(saveImageAction, SIGNAL(triggered()), this, SLOT(saveImage()));

  QAction *deleteImageAction = new QAction(tr("Delete image"), this);
  deleteImageAction->setIcon(QIcon(":/Images/famfamfam_silk_icons/picture_delete.png"));
  connect(deleteImageAction, SIGNAL(triggered()), this, SLOT(deleteImage()));

  setContextMenuPolicy(Qt::ActionsContextMenu);
  addAction(saveImageAction);
  addAction(deleteImageAction);
}

MorphologyImageWidget::~MorphologyImageWidget() {
}

void MorphologyImageWidget::resizeEvent(QResizeEvent *event) {
  int width = this->width(); //parentWidget()->height();
  int height = (width * pixmap_.height()) / pixmap_.width();
  setPixmap(pixmap_.scaled(width, height));

  setFixedHeight(height);

  QLabel::resizeEvent(event);
}


// private slot
void MorphologyImageWidget::saveImage() {
  ImageSaver::savePixmap(pixmap_, this);
}

// private slot
void MorphologyImageWidget::deleteImage() {
  setVisible(false);
  emit deleted();
}

}