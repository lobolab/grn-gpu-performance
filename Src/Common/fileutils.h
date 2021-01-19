// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#pragma once

#include <QFile>

namespace LoboLab {

namespace FileUtils {
// QFile::copy also copy the attributes, which is not right
// Solution: let's make a completely new file
static void copyRaw(const QString &srcFileN, const QString &newFileN) {
  
  if (srcFileN != newFileN) {
    QFile::remove(newFileN);

    QFile newFile(newFileN);
    QFile srcFile(srcFileN);
    newFile.open(QFile::WriteOnly);
    srcFile.open(QFile::ReadOnly);

    uint dataLength = 16000;
    char *data = new char[dataLength];
    while (!srcFile.atEnd()) {
      qint64 len = srcFile.read(data, dataLength);
      newFile.write(data, len);
    }

    delete []data;
    newFile.close();
    srcFile.close();
  }
}
};

}