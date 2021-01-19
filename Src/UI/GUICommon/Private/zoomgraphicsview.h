// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#pragma once

#include <QGraphicsView>

namespace LoboLab {

class ZoomGraphicsView : public QGraphicsView {
  Q_OBJECT

   public:
    ZoomGraphicsView(QWidget *parent = NULL, bool keepFitted = true,
                     bool reversed = true);
    virtual ~ZoomGraphicsView();

    inline QGraphicsScene *getScene() { return scene_; }

    void fitView();

   public slots:
    void saveImage();

   protected:
    virtual void wheelEvent(QWheelEvent *e);
    virtual void keyPressEvent(QKeyEvent *e);
    virtual void keyReleaseEvent(QKeyEvent *e);
    virtual void focusOutEvent(QFocusEvent *e);
    virtual void resizeEvent(QResizeEvent * event);

    QGraphicsScene *scene_;
    QAction *saveImageAction_;

   private:

    bool keepFitted_;
    bool reversed_;
    bool zoomMode_;
    bool fittingView_;
};

} // namespace LoboLab
