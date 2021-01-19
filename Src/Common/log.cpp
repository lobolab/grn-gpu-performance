// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#include "log.h"

#include <QFile>
#include <QTextStream>
#include <QDate>

namespace LoboLab {

Log *Log::singleton_ = NULL;

Log::Log() : file_(NULL), stream_(NULL), lockingStream_(NULL) {}

Log::~Log() {
  delete lockingStream_;
  delete stream_;
  delete file_;
}

void Log::setLogName(const QString &fileName) {
  if (!singleton_)
    singleton_ = new Log();
  else {
    delete singleton_->lockingStream_;
    delete singleton_->stream_;
    delete singleton_->file_;
  }

  singleton_->file_ = new QFile(fileName);
}

void Log::openLog() {
  singleton_->file_->open(QIODevice::Append | QIODevice::WriteOnly |
                        QIODevice::Text);

  Q_ASSERT(singleton_->file_->isOpen());

  singleton_->stream_ = new QTextStream(singleton_->file_);
  singleton_->lockingStream_ = new LockingTextStream(*singleton_->stream_);

  Log::write() << "###################### LOG START ########################"
               << endl;
}

void Log::closeLog() {
  if (singleton_) {
    singleton_->file_->close();
    delete singleton_;
    singleton_ = NULL;
  }
}

Log::LockingTextStream &Log::write() {
  if (!singleton_) {
    setLogName("out.log");
    singleton_->openLog();
  } else if (!singleton_->file_->isOpen())
    singleton_->openLog();

  return *singleton_->lockingStream_ << QDateTime::currentDateTimeUtc()
         .toString("yyyy-MM-ddThh:mm:ss.zzzZ : ");
}

}