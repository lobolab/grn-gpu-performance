// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#include "colorsconfig.h"

namespace LoboLab {

const int ColorsConfig::nColors = 19;

// pure colors (http://www.quackit.com/html/codes/color/color_names.cfm)
const QRgb ColorsConfig::colorsMorp[nColors] = {
  0xFF0000, //red
  0x00FF00, // green
  0x0000FF, // blue
  0xFF00FF, // fuchsia
  0xFFFF00, // yellow
  0x00FFFF, // cyan
  0xC0C0C0, // silver
  0x7425B1, // purple dark
  0xFF9233, // orange
  0xAD2323, // dark red
  0xFED4D9, // pink
  0x5F6142, // artichoke
  0x43B3AE, // verdigris
  0xFEDA83, // saffron
  0x814A19, // brown
  0xC0D8F2, // bashful blue
  0xC3C982, // celery
  0xE88F2B, // more mustard
  0x868A35, // old olive
};

// crayonbow (http://colrd.com/create/palette/?id=18733)
const QRgb ColorsConfig::colorsCells[nColors] = {
  0xF6300A, //red
  0x00C000, // green
  0x0A62DA, // blue
  0xFF00FF, // fuchsia
  0xFFFF00, // yellow
  0x00FFFF, // cyan
  0xC0C0C0, // silver
  0x7425B1, // purple dark
  0xFF9233, // orange
  0xAD2323, // dark red
  0xFED4D9, // pink
  0x5F6142, // artichoke
  0x43B3AE, // verdigris
  0xFEDA83, // saffron
  0x814A19, // brown
  0xC0D8F2, // bashful blue
  0xC3C982, // celery
  0xE88F2B, // more mustard
  0x868A35, // old olive
};

const QRgb ColorsConfig::background = 0x00000000;

}