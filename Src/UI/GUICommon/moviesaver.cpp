// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#include "moviesaver.h"
#include "Common/mathalgo.h"

#include <QDir>
#include <QFileDialog>
#include <QPixmap>
#include <QWidget>
//#include <QSvgGenerator>
#include <QPainter>
#include <QMessageBox>
#include <QSettings>
#include <QProcess>

namespace LoboLab {

MovieSaver::MovieSaver(QWidget *parent) 
  : parent_(parent), 
    recordingMovie_(false) {
  tempDir_ = "temp_movie_frames";

}

MovieSaver::~MovieSaver() {
  cancelMovie();
}

void MovieSaver::startMovie(QWidget *widgetToRecord) {
  if (isRecording())
    cancelMovie();

  widgetToRecord_ = widgetToRecord;
  recordingMovie_ = true;

  QDir dir;
  dir.mkdir(tempDir_);
  iFrame_ = 0;
  saveFrame();
}


void MovieSaver::saveFrame() {
  if (recordingMovie_) {
    //QPixmap pixmap = QPixmap::grabWidget(widgetToRecord_);
    QPixmap pixmap = QPixmap::grabWindow(widgetToRecord_->winId());
    pixmap.save(tempDir_ + QString("/frame%1.png").arg(iFrame_, 5, 10,
      QChar('0')));
    iFrame_++;
  }
}

void MovieSaver::endMovie(QString fileName) {
  if (recordingMovie_) {
    if (fileName.isEmpty()) {
      QSettings settings;
      settings.beginGroup("MovieSaver");
      QString lastDir = settings.value("lastDir").toString();

      fileName = QFileDialog::getSaveFileName(parent_, "Save movie", lastDir, 
                                              "Movies (*.mp4)");

      if (!fileName.isEmpty()) {
        QString newDir = QFileInfo(fileName).absolutePath();
        if (newDir != lastDir)
          settings.setValue("lastDir", newDir);

        QFile::remove(fileName);
      }
    }

    if (!fileName.isEmpty()) {
      QProcess process(parent_);

      process.start("ffmpeg -i " + tempDir_ + "/frame%05d.png -c:v libx264 "
                    "-pix_fmt yuv444p -crf 0 " + fileName);
      process.waitForFinished();
    }
      
    removeTempFiles();
    widgetToRecord_ = NULL;
    iFrame_ = 0;
    recordingMovie_ = false;
  }
}

void MovieSaver::cancelMovie() {
  if (recordingMovie_) {
    widgetToRecord_ = NULL;
    iFrame_ = 0;
    recordingMovie_ = false;
    removeTempFiles();
  }
}

void MovieSaver::removeTempFiles() {
  QDir dir(QDir::currentPath() + "/" + tempDir_);
  if (dir.exists(QDir::currentPath() + "/" + tempDir_)) {
    Q_FOREACH(QFileInfo info, dir.entryInfoList(QDir::NoDotAndDotDot | QDir::System | QDir::Hidden  | QDir::AllDirs | QDir::Files, QDir::DirsFirst)) {
        Q_ASSERT(!info.isDir());
        QFile::remove(info.absoluteFilePath());
    }
    dir.rmdir(QDir::currentPath() + "/" + tempDir_);
  }
}

}
