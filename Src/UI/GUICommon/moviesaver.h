// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#pragma once

#include <QWidget>

namespace LoboLab {

class MovieSaver {
 public:
  MovieSaver(QWidget *parent);
  ~MovieSaver();

  void startMovie(QWidget *widgetToRecord);
  void saveFrame();
  void endMovie(QString fileName = QString());
  void cancelMovie();

  inline bool isRecording() { return recordingMovie_; }

 private:
  void removeTempFiles();

  QWidget *parent_;
  QWidget *widgetToRecord_;
  int iFrame_;
  bool recordingMovie_;
  QString tempDir_;

};

} // namespace LoboLab
