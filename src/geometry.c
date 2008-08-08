/** \file geometry.h           Functions for calculating distances and angles.
 *
 *             Part of MDT, Copyright(c) 1989-2008 Andrej Sali
 */

#include <glib.h>
#include <math.h>
#include "modeller.h"
#include "geometry.h"
#include "mdt_index.h"

/** Return the distance between two coordinates */
static float dist1(float x1, float y1, float z1, float x2, float y2, float z2)
{
  float xd, yd, zd;
  xd = x1 - x2;
  yd = y1 - y2;
  zd = z1 - z2;
  return sqrt(xd * xd + yd * yd + zd * zd);
}

/** Return the bin index for the distance between two specified atoms in the
    same protein. */
int idist0(int ia1, int ia1p, const struct mod_structure *struc,
           const struct mod_mdt_libfeature *feat)
{
  if (ia1 >= 0 && ia1p >= 0) {
    float d, *x, *y, *z;
    x = mod_float1_pt(&struc->cd.x);
    y = mod_float1_pt(&struc->cd.y);
    z = mod_float1_pt(&struc->cd.z);
    d = dist1(x[ia1], y[ia1], z[ia1], x[ia1p], y[ia1p], z[ia1p]);
    return iclsbin(d, feat);
  } else {
    return feat->nbins;
  }
}
