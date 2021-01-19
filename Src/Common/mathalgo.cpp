// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#include "mathalgo.h"

namespace LoboLab {

// rotate the point (pX,pY) aroung the point (oX, oY) for ang radians.
void MathAlgo::rotate(double &pX, double &pY,
                      const double oX, const double oY, const double ang) {
  const double s = sin(ang);
  const double c = cos(ang);

  // translate point back to origin:
  pX -= oX;
  pY -= oY;

  // rotate point
  const double tX = pX * c + pY * s;
  const double tY = -pX * s + pY * c;

  // translate point back:
  pX = tX + oX;
  pY = tY + oY;
}

// rotate the point (pX,pY) aroung the point (oX, oY) for ang radians.
void MathAlgo::rotate(int &pX, int &pY, const int oX, const int oY,
                      const double ang) {
  const double s = sin(ang);
  const double c = cos(ang);

  // translate point back to origin:
  pX -= oX;
  pY -= oY;

  // rotate point
  const double tX = pX * c + pY * s;
  const double tY = -pX * s + pY * c;

  // translate point back:
  pX = tX + oX;
  pY = tY + oY;
}

// minimum distance between point p and segment vw
double MathAlgo::distSegment(const QPointF &p, const QPointF &v,
                             const QPointF &w) {
  const double l2 = distSquared(v, w);
  if (l2 == 0.0) // a == b)
    return dist(p, v);

  // Consider the line extending the segment, parameterized as
  // v + t (w - v).
  // We find projection of point p onto the line.
  // It falls where t = [(p-v) . (w-v)] / |w-v|^2
  const double t = dot(p - v, w - v) / l2;
  if (t < 0.0)		 // Beyond the 'v' end of the segment
    return dist(p, v);
  else if (t > 1.0) // Beyond the 'w' end of the segment
    return dist(p, w);

  const QPointF projection = v + t * (w - v);
  return dist(p, projection);
}

// See http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html
bool MathAlgo::insideArea(const int pX, const int pY, const QVector<QPointF> &area) {
  int i, j;
  bool c = false;
  int nPoints = area.size();
  for (i = 0, j = nPoints-1; i < nPoints; j = i++) {
    if ( ((area[i].y() > pY) != (area[j].y() > pY)) &&
         (pX < (area[j].x() - area[i].x()) *
          (pY - area[i].y()) / (area[j].y() -
                                area[i].y()) + area[i].x()) ) {
      c = !c;
    }
  }

  return c;
}


} // namespace LoboLab
