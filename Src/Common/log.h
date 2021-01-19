// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#pragma once

#include <QTextStream>
#include <QFile>
#include <QMutex>

namespace LoboLab {

class Log {
 public:
  ~Log();

  class LockingTextStream {
   public:
    LockingTextStream(QTextStream &s)
      : stream_(s) {
    }

    template<typename T> inline
    LockingTextStream& operator<<(const T &t) {
      mutex_.lock();
      stream_ << t;
      mutex_.unlock();
      return *this;
    }

    QMutex mutex_;
    QTextStream &stream_;
  };

  static void setLogName(const QString &fileName);
  static void closeLog();
  static LockingTextStream &write();

 private:
  Log();
  Q_DISABLE_COPY(Log)

  void openLog();

  static Log *singleton_;
  QFile *file_;
  QTextStream *stream_;
  LockingTextStream *lockingStream_;
};

} // namespace LoboLab