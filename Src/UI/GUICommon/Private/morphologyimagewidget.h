// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#pragma once

#include <QLabel>

namespace LoboLab {

class MorphologyImage;

class MorphologyImageWidget : public QLabel {
  Q_OBJECT

  public:
    MorphologyImageWidget(MorphologyImage *morphologyImage, QWidget *parent = NULL);
    virtual ~MorphologyImageWidget();

    inline const QPixmap &getPixmap() const { return pixmap_; }

  signals:
    void deleted();

  protected:
    virtual void resizeEvent(QResizeEvent *event);

    private slots:
    void saveImage();
    void deleteImage();

  private:
    MorphologyImage *morphologyImage_;
    QPixmap pixmap_;
  };

} // namespace LoboLab
