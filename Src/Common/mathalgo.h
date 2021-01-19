// Copyright (c) Lobo Lab (lobolab.umbc.edu)
// All rights reserved.

#pragma once

#if !defined(_USE_MATH_DEFINES)
#define _USE_MATH_DEFINES 1
#endif

//#include <math.h>
#include <QtMath>

// Very fast math operations for amd machines
#ifdef USE_AMDLIBM
#define REPLACE_WITH_AMDLIBM
#include "amdlibm.h"
#endif

#include <Eigen/Core> // here to be sure it is defined after AMDLIBM

namespace Eigen {
  typedef Array<bool, Eigen::Dynamic, Eigen::Dynamic> MatrixXb;
}

#include <random>

#include <QVector>
#include <QPointF>

namespace LoboLab {

namespace MathAlgo {


inline bool checkBit(int var, int pos) {
  return var & (1<<pos);
}

// From: http://martin.ankerl.com/2012/01/25/optimized-approximative-pow-in-c-and-cpp/
inline double fastPow(double a, double b) {
  union {
    double d;
    int x[2];
  } u = { a };
  u.x[1] = (int)(b * (u.x[1] - 1072632447) + 1072632447);
  u.x[0] = 0;
  return u.d;
}

inline int round(double n) {
  return n-floor(n)>=0.5? ceil(n) : floor(n);
}

// d significant digits. See http://stackoverflow.com/questions/13094224/a-c-routine-to-round-a-float-to-n-significant-digits
inline double ceilS(double n, int d) {
  if (n == 0.0)
    return 0.0;
  else {
    double factor = pow(10.0, d - ceil(log10(fabs(n))));
    return ceil(n * factor) / factor;   
  }
}

inline double dist(double aX, double aY, double bX, double bY) {
  return sqrt(pow(bX-aX, 2) + pow(bY-aY, 2));
}

inline double dist(const QPointF &a, const QPointF &b) {
  return sqrt(pow(b.x()-a.x(), 2) + pow(b.y()-a.y(), 2));
}

inline double dist(const QPoint &a, const QPoint &b) {
  return sqrt(pow((double)b.x()-a.x(), 2) + pow((double)b.y()-a.y(), 2));
}

inline double dist(const double a, const double b) {
  return qAbs(b-a);
}

// square of the distance between a and b (faster than calculating the dist)
inline double distSquared(const QPointF &a, const QPointF &b) {
  return pow(b.x()-a.x(), 2) + pow(b.y()-a.y(), 2);
}

// rotate the point (pX,pY) aroung the point (oX, oY) for ang radians.
void rotate(double &pX, double &pY,
            const double oX, const double oY, const double ang);

// rotate the point (pX,pY) aroung the point (oX, oY) for ang radians.
void rotate(int &pX, int &pY, const int oX, const int oY,
            const double ang);

inline double length(const QPointF &v) {
  return sqrt(pow(v.x(), 2) + pow(v.y(), 2));
}

inline double length(double vX, double vY) {
  return sqrt(pow(vX, 2) + pow(vY, 2));
}

inline double dot(const QPointF &a, const QPointF &b) {
  return a.x() * b.x() + a.y() * b.y();
}

// minimum distance between point p and segment vw
double distSegment(const QPointF &p, const QPointF &v,
                   const QPointF &w);

inline double deg2rad(double ang) {
  return ang*M_PI/180;
}

inline double rad2deg(double ang) {
  return ang*180/M_PI;
}

inline double diffAngle(double ang1, double ang2) {
  double d = atan2(sin((ang2-ang1)), cos((ang2-ang1)));
  if (d < 0)
    d += 2*M_PI;

  return d;
}

inline double midAngle(double ang1, double ang2) {
  double d = MathAlgo::diffAngle(ang1, ang2);
  double m = ang1 + d / 2;
  return m;
}

inline double invAngle(double ang) {
  if (ang>=M_PI)
    return ang - M_PI;
  else
    return ang + M_PI;
}

// Positive modulo, b must be positive
inline int modPos(int a, int b) {
  int ret = a % b;
  if (ret < 0)
    ret += b;
  return ret;
}

inline double modDouble(double a, int b) {
  return a - floor(a/b) * b;
}

inline double modAng(double ang) {
  return ang - floor(ang/(2*M_PI)) * 2*M_PI;
}

inline double minDiffAngle(double ang1, double ang2) {
  double d = modAng(ang1-ang2);
  if (d > M_PI)
    d = 2*M_PI - d;

  return d;
}

inline QPointF vector(double dist, double ang) {
  return QPointF(dist*cos(ang), dist*sin(ang));
}

inline QPointF dirOfVector(const QPointF &vect) {
  double vectLen = length(vect);
  if (vectLen>0.0001)
    return vect/vectLen;
  else
    return QPointF();
}

inline QPointF dirOfPoints(const QPointF &from, const QPointF &to) {
  return (to-from)/dist(from, to);
}

inline void dirOfAng(double &dirX, double &dirY, double ang) {
  dirX = cos(ang);
  dirY = sin(ang);
}

inline QPointF dirOfAng(const double ang) {
  return QPointF(cos(ang), sin(ang));
}

inline double angOfDir(const double dirX, const double dirY) {
  double ang = atan2(dirY, dirX);
  if (ang < 0)
    ang += 2*M_PI;

  return ang;
}

inline double angOfDir(const QPointF &dir) {
  double ang = atan2(dir.y(), dir.x());
  if (ang < 0)
    ang += 2*M_PI;

  return ang;
}

inline QPointF midDir(const QPointF &a, const QPointF &b) {
  return MathAlgo::dirOfAng(MathAlgo::midAngle(
                              MathAlgo::angOfDir(a), MathAlgo::angOfDir(b)));
}

template <class T>
inline const T& max(const T& a, const T& b) {
  return (a<b)?b:a;
}

template <class T>
inline const T& min(const T& a, const T& b) {
  return (!(b<a))?a:b;
}

// See http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html
bool insideArea(const int pX, const int pY, const QVector<QPointF> &area);

// Returns a different number each milisecond 
inline unsigned int randSeed() {
  static std::random_device rdev;
  return rdev();

  //return  1000 * QDateTime::currentDateTimeUtc().toTime_t() +
  //        (uint) QTime::currentTime().msec();
}

static std::mt19937_64 rand_gen(randSeed());

inline bool randBool() {
  return (rand_gen() % 2) == 0;
}

// rand integer between 0 and n-1 (suffer small errors)
inline int randInt(int n) {
  return rand_gen() % n;
}

// rand integer between min and max (inclusive)
inline int randInt(int min, int max) {
  return min + (rand_gen() % (max - min + 1));
}

inline int rand100() {
  static std::uniform_int_distribution<int> dist(0, 99);
  return dist(rand_gen);
}

inline int rand1000() {
  static std::uniform_int_distribution<int> dist(0, 999);
  return dist(rand_gen);
}

inline double rand01() {
  static std::uniform_real_distribution<double> dist(0.0, 1.0);
  return dist(rand_gen);
}

inline double rand15() {
  static std::uniform_real_distribution<double> dist(1.0, 5.0);
  return dist(rand_gen);
}

// Shuffle algorimth: Durstenfeld, Richard (July 1964).
// "Algorithm 235: Random permutation". Communications of the ACM 7 (7):
// 420. doi:10.1145/364520.364540.
template <class T>
void shuffle(int n, T *vec) {
  for (int i = n-1; i > 0; --i) {
    int j = rand_gen() % (i+1);
    T t = vec[j];
    vec[j] = vec[i];
    vec[i] = t;
  }
}

template <class T>
void shuffle(QList<T> &list) {
  int n = list.size();
  for (int i = n-1; i > 0; --i) {
    int j = rand_gen() % (i+1);
    list.swap(i,j);
  }
}

} // namespace MathAlgo

} // namespace LoboLab


