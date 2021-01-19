// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#pragma once

#include "Common/mathalgo.h"
#include <QColor>

namespace LoboLab {

namespace GuiConfig {

  inline double paramFactor() { return 1; }
  inline double regionDistFactor() { return 1; }
  inline double organDistFactor() { return 1; }

  inline double eyeSize() { return 20; }
  inline double brainSize() { return 40; }
  inline double pharynxSize() { return 60;}
  inline double VNCWidth() { return 8; }

  inline double organNodeRadius() { return 15; }
  inline double regionNodeRadius() { return 20; }
  inline double nodeTextPadding() { return 5; }
  inline double regionLeverRadius() { return 5; }
  inline double linkLeverRadius() { return 5; }
  inline double organLeverRadius() { return 3.5; }
  inline double MGLineOrganEndRadius() { return 3.5; }
  inline double graphRegionsLinkWidth() { return 1.5; }
  inline double graphOrgansLinkWidth() { return 1; }

  inline const QColor &regionLeverColor() {
    static QColor color(255, 0, 0);
    return color;
  }

  inline const QColor &linkLeverColor() {
    static QColor color(0, 255, 0);
    return color;
  }

  inline const QColor &organLeverColor() {
    static QColor color(0, 0, 255);
    return color;
  }

  inline const QColor &graphLinkColor() {
    static QColor color(Qt::darkRed);
    return color;
  }

  inline const QColor &MGLineOrganEndColor() {
    static QColor color(Qt::darkYellow);
    return color;
  }

  inline const QColor &cropAreaPointColor() {
    static QColor color(255, 255, 0);
    return color;
  }

  inline const QColor &radiationColor() {
    static QColor color(255, 254, 2, 150);
    return color;
  }

  inline double defaultRegionsLinkDist() { return 200; }
  inline double defaultRegionsLinkAng() { return 3*M_PI/2; }
  inline double defaultOrganRot() { return 0; }
  inline double defaultOrganDist() { return 20; }
  inline double defaultOrganAng() { return 0; }
  inline double defaultLineOrganEndDist() { return 20; }
  inline double defaultLineOrganEndAng() { return M_PI_2; }

  inline double phenRoundness() { return 0.25; }

  inline const QPointF &actionNodeYDisp() {
    static QPointF vector(0, 75);
    return vector;
  }
  inline const QPointF &actionNodeXDisp() {
    static QPointF vector(50, 0);
    return vector;
  }
  inline double actionNodeRadius() { return 22; }
  inline double cropAreaPointRadius() { return 5; }
} // namespace GuiConfig

} // namespace LoboLab
