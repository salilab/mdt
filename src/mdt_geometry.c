/** \file mdt_geometry.c  Functions to calculate geometrical properties
 *
 *             Part of MDT, Copyright(c) 1989-2007 Andrej Sali
 */

#include <math.h>
#include <glib.h>

/** Return the distance between two coordinates */
float dist1(float x1, float y1, float z1, float x2, float y2, float z2)
{
  float xd, yd, zd;
  xd = x1 - x2;
  yd = y1 - y2;
  zd = z1 - z2;
  return sqrt(xd * xd + yd * yd + zd * zd);
}

/** Return the angle between three coordinates */
float angle1(float x1, float y1, float z1, float x2, float y2, float z2,
             float x3, float y3, float z3)
{
  static const float tiny = 1.0e-15;
  float d1, d2, v1x, v1y, v1z, v2x, v2y, v2z, scalprod, sizeprod, div;
  v1x = x1 - x2;
  v1y = y1 - y2;
  v1z = z1 - z2;
  v2x = x3 - x2;
  v2y = y3 - y2;
  v2z = z3 - z2;
  d1 = sqrt(v1x * v1x + v1y * v1y + v1z * v1z);
  d2 = sqrt(v2x * v2x + v2y * v2y + v2z * v2z);
  scalprod = v1x * v2x + v1y * v2y + v1z * v2z;
  sizeprod = d1 * d2;
  div = (sizeprod > tiny ? scalprod / sizeprod : 0.0);
  div = CLAMP(div, -1.0, 1.0);
  return acos(div) * 180.0 / G_PI;
}

/** Return the dihedral angle between four coordinates */
float dihedral1(float x1, float y1, float z1, float x2, float y2,
                float z2, float x3, float y3, float z3, float x4,
                float y4, float z4)
{
  double rt[4][3], l1[3], l2[3], l3[3], xt1[3], xt2[3], leng1, leng2, dot1,
      ang, sign, norm;

  rt[0][0] = x1;
  rt[0][1] = y1;
  rt[0][2] = z1;
  rt[1][0] = x2;
  rt[1][1] = y2;
  rt[1][2] = z2;
  rt[2][0] = x3;
  rt[2][1] = y3;
  rt[2][2] = z3;
  rt[3][0] = x4;
  rt[3][1] = y4;
  rt[3][2] = z4;

  l1[0] = rt[1][0] - rt[0][0];
  l1[1] = rt[1][1] - rt[0][1];
  l1[2] = rt[1][2] - rt[0][2];
  l2[0] = rt[2][0] - rt[1][0];
  l2[1] = rt[2][1] - rt[1][1];
  l2[2] = rt[2][2] - rt[1][2];
  l3[0] = rt[3][0] - rt[2][0];
  l3[1] = rt[3][1] - rt[2][1];
  l3[2] = rt[3][2] - rt[2][2];

  xt1[0] = l2[1] * l1[2] - l2[2] * l1[1];
  xt1[1] = l2[2] * l1[0] - l2[0] * l1[2];
  xt1[2] = l2[0] * l1[1] - l2[1] * l1[0];
  xt2[0] = l3[1] * l2[2] - l3[2] * l2[1];
  xt2[1] = l3[2] * l2[0] - l3[0] * l2[2];
  xt2[2] = l3[0] * l2[1] - l3[1] * l2[0];
  leng1 = xt1[0] * xt1[0] + xt1[1] * xt1[1] + xt1[2] * xt1[2];
  leng2 = xt2[0] * xt2[0] + xt2[1] * xt2[1] + xt2[2] * xt2[2];
  dot1 = xt1[0] * xt2[0] + xt1[1] * xt2[1] + xt1[2] * xt2[2];
  norm = dot1 / sqrt(leng1 * leng2);
  norm = CLAMP(norm, -1.0, 1.0);
  ang = acos(norm);
  sign = xt1[0] * l3[0] + xt1[1] * l3[1] + xt1[2] * l3[2];
  if (sign < 0.0) {
    ang = -ang;
  }
  return -ang * 180.0 / G_PI;
}
