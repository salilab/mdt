/** \file mdt_index.c  Functions to calculate MDT indices.
 *
 *             Part of MDT, Copyright(c) 1989-2007 Andrej Sali
 */

#include <string.h>
#include <math.h>
#include "modeller.h"
#include "util.h"
#include "mdt_index.h"
#include "mdt_property.h"
#include "mdt_hydrogen_bonds.h"
#include "mdt_stereo.h"
#include "mdt_tuples.h"

static int ires_get(int ires, int nres, int igaptyp, const int *irestyp,
                    int ndimen)
{
  if (ires < 1 || ires > nres) {
    return igaptyp;
  } else if (irestyp[ires - 1] >= 21 || irestyp[ires - 1] <= 0) {
    mod_logwarning("irestab", "ires, irestp: %d %d", ires, irestyp[ires - 1]);
    return ndimen;
  } else {
    return irestyp[ires - 1];
  }
}

static int irestab(const struct f_int2_array *ialn, int naln, int iseq,
                   const int *irestyp, int nres, int ip, int delta,
                   gboolean delta_ali, int ndimen, int igaptyp)
{
  if (delta_ali) {
    int ipos = ip + delta;
    if (ipos < 0 || ipos >= naln) {
      return igaptyp;
    } else {
      int ires = f_int2_get(ialn, ipos, iseq);
      return ires_get(ires, nres, igaptyp, irestyp, ndimen);
    }
  } else {
    int ires = f_int2_get(ialn, ip, iseq);
    if (ires < 1 || ires > nres) {
      return igaptyp;
    } else {
      return ires_get(ires + delta, nres, igaptyp, irestyp, ndimen);
    }
  }
}

/** Return the bin index itab[ir], or the undefined bin index if anything is
    out of range */
static int itable(const int *itab, int nr, int ir, int ndim)
{
  if (ir >= 0 && ir < nr && itab[ir] >= 1 && itab[ir] <= ndim) {
    return itab[ir];
  } else {
    return ndim;
  }
}

/** Ensure that a given bin index is in range for the feature; return in
    the undefined bin if not. */
static int index_inrange(int index, const struct mdt_libfeature *feat)
{
  return (index >= 1 && index < feat->nbins) ? index : feat->nbins;
}

/** Convert a raw number to the corresponding MDT bin index */
int iclsbin(float x, const struct mdt_library *mlib, int ifi, int nrang)
{
  int i;
  const struct mdt_libfeature *feat = &mlib->base.features[ifi - 1];
  const struct mdt_bin *bin = feat->bins;
  for (i = 0; i < nrang - 1; i++, bin++) {
    if (x >= bin->rang1 && x <= bin->rang2) {
      return i + 1;
    }
  }
  bin = &feat->bins[0];
  mod_logwarning("iclsbin", "Undefined value; X,x1,x2,n,bin: %f %f %f %d",
                 x, bin->rang1, bin->rang2, nrang);
  return nrang;
}

/** Return the distance between two coordinates */
static float dist1(float x1, float y1, float z1, float x2, float y2, float z2)
{
  float xd, yd, zd;
  xd = x1 - x2;
  yd = y1 - y2;
  zd = z1 - z2;
  return sqrt(xd * xd + yd * yd + zd * zd);
}

/** Return the angle between three coordinates */
static float angle1(float x1, float y1, float z1, float x2, float y2, float z2,
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

/** Return the dihedral angle between four coordinates.
    outrange is set to TRUE if the angle cannot be reliably calculated. */
static float dihedral1(float x1, float y1, float z1, float x2, float y2,
                       float z2, float x3, float y3, float z3, float x4,
                       float y4, float z4, gboolean *outrange)
{
  double rt[4][3], l1[3], l2[3], l3[3], xt1[3], xt2[3], leng1, leng2, dot1,
      ang, sign, norm, lengprod;

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
  lengprod = leng1 * leng2;
  if (lengprod < 1.0e-4) {
    *outrange = TRUE;
    return 0.;
  } else {
    *outrange = FALSE;
    norm = dot1 / sqrt(lengprod);
  }
  norm = CLAMP(norm, -1.0, 1.0);
  ang = acos(norm);
  sign = xt1[0] * l3[0] + xt1[1] * l3[1] + xt1[2] * l3[2];
  if (sign < 0.0) {
    ang = -ang;
  }
  return -ang * 180.0 / G_PI;
}


/** Return the bin index for the distance between two specified atoms in the
    same protein. */
static int idist0(int ia1, int ia1p, const struct structure *struc,
                  const struct mdt_library *mlib, int ifi, int nrang)
{
  if (ia1 >= 0 && ia1p >= 0) {
    float d, *x, *y, *z;
    x = f_float1_pt(&struc->cd.x);
    y = f_float1_pt(&struc->cd.y);
    z = f_float1_pt(&struc->cd.z);
    d = dist1(x[ia1], y[ia1], z[ia1], x[ia1p], y[ia1p], z[ia1p]);
    return iclsbin(d, mlib, ifi, nrang);
  } else {
    return nrang;
  }
}

/** Return the bin index for the angle between three specified atoms in the
    same protein. */
static int iangle0(int ia1, int ia2, int ia3, const struct structure *struc,
                   const struct mdt_library *mlib, int ifi, int nrang)
{
  if (ia1 >= 0 && ia2 >= 0 && ia3 >= 0) {
    float d, *x, *y, *z;
    x = f_float1_pt(&struc->cd.x);
    y = f_float1_pt(&struc->cd.y);
    z = f_float1_pt(&struc->cd.z);
    d = angle1(x[ia1], y[ia1], z[ia1], x[ia2], y[ia2], z[ia2], x[ia3],
               y[ia3], z[ia3]);
    return iclsbin(d, mlib, ifi, nrang);
  } else {
    return nrang;
  }
}

/** Return the bin index for the dihedral angle between four specified atoms
    in the same protein. */
static int idihedral0(int ia1, int ia2, int ia3, int ia4,
                      const struct structure *struc,
                      const struct mdt_library *mlib, int ifi, int nrang)
{
  if (ia1 >= 0 && ia2 >= 0 && ia3 >= 0 && ia4 >= 0) {
    gboolean outrange;
    float d, *x, *y, *z;
    x = f_float1_pt(&struc->cd.x);
    y = f_float1_pt(&struc->cd.y);
    z = f_float1_pt(&struc->cd.z);
    d = dihedral1(x[ia1], y[ia1], z[ia1], x[ia2], y[ia2], z[ia2], x[ia3],
                  y[ia3], z[ia3], x[ia4], y[ia4], z[ia4], &outrange);
    if (outrange) {
      return nrang;
    } else {
      return iclsbin(d, mlib, ifi, nrang);
    }
  } else {
    return nrang;
  }
}

/** Get the index into the MDT for the given alignment feature */
int my_mdt_index(int ifi, const struct alignment *aln, int is1, int ip1,
                 int is2, int ir1, int ir2, int ir1p, int ir2p, int ia1,
                 int ia1p, const struct mdt_library *mlib, int ip2,
                 int ibnd1, int ibnd1p, int is3, int ir3, int ir3p,
                 const struct libraries *libs,
                 const struct energy_data *edat,
                 struct mdt_properties *prop, GError **err)
{
  int ret, ierr = 0;
  const int *binprop;
  int iresol, irad;
  float fprop;
  struct structure *struc1, *struc2;
  struct sequence *seq1, *seq2;
  const struct mdt_bond *bond;
  const struct mdt_tuple *tup, *tup2;
  struct mdt_libfeature *feat = &mlib->base.features[ifi - 1];
  struc1 = alignment_structure_get(aln, is1);
  struc2 = alignment_structure_get(aln, is2);
  seq1 = alignment_sequence_get(aln, is1);
  seq2 = alignment_sequence_get(aln, is2);
  switch (ifi) {
  case 35:
    iresol = property_iresol(aln, is1, prop, mlib, ifi, feat);
    return index_inrange(iresol, feat);
  case 38:
    iresol = property_iresol(aln, is2, prop, mlib, ifi, feat);
    return index_inrange(iresol, feat);
  case 66:
    return irestab(&aln->ialn, aln->naln, is1, f_int1_pt(&seq1->irestyp),
                   seq1->nres, ip1, mlib->deltai, mlib->deltai_ali,
                   feat->nbins, libs->igaptyp);
  case 67:
    return irestab(&aln->ialn, aln->naln, is2, f_int1_pt(&seq2->irestyp),
                   seq2->nres, ip1, mlib->deltai, mlib->deltai_ali,
                   feat->nbins, libs->igaptyp);
  case 77:
    return irestab(&aln->ialn, aln->naln, is1, f_int1_pt(&seq1->irestyp),
                   seq1->nres, ip1, mlib->deltaj, mlib->deltaj_ali,
                   feat->nbins, libs->igaptyp);
  case 78:
    return irestab(&aln->ialn, aln->naln, is2, f_int1_pt(&seq2->irestyp),
                   seq2->nres, ip1, mlib->deltaj, mlib->deltaj_ali,
                   feat->nbins, libs->igaptyp);
  case 79:
    binprop = property_iatta(aln, is1, prop, mlib, ifi, libs, err);
    if (!binprop) {
      return 0.;
    }
    return itable(binprop, struc1->cd.natm, ia1, feat->nbins);
  case 80:
    return itable(property_iatmacc(aln, is1, prop, mlib, ifi, feat),
                  struc1->cd.natm, ia1, feat->nbins);
  case 81:
    binprop = property_ifatmacc(aln, is1, prop, mlib, ifi, feat, libs, err);
    if (!binprop) {
      return 0;
    }
    return itable(binprop, struc1->cd.natm, ia1, feat->nbins);
  case 82:
  case 103:
    return idist0(ia1, ia1p, struc1, mlib, ifi, feat->nbins);
  case 83:
    binprop = property_iatta(aln, is1, prop, mlib, ifi, libs, err);
    if (!binprop) {
      return 0;
    }
    return itable(binprop, struc1->cd.natm, ia1p, feat->nbins);
  case 84:
    binprop = property_hb_iatta(aln, is1, prop, mlib, ifi, libs, err);
    if (!binprop) {
      return 0;
    }
    return numb_hda(ia1, binprop, &struc1->cd, mlib->hbond, mlib->hbond_cutoff,
                    0, feat->nbins);
  case 85:
    binprop = property_hb_iatta(aln, is1, prop, mlib, ifi, libs, err);
    if (!binprop) {
      return 0;
    }
    return numb_hda(ia1, binprop, &struc1->cd, mlib->hbond, mlib->hbond_cutoff,
                    1, feat->nbins);
  case 86:
    if (!property_hbpot(aln, is1, prop, mlib, ifi, libs, &fprop, err)) {
      return 0;
    }
    return iclsbin(fprop, mlib, ifi, feat->nbins);
  case 87:
    binprop = property_hb_iatta(aln, is1, prop, mlib, ifi, libs, err);
    if (!binprop) {
      return 0;
    }
    return numb_hda(ia1, binprop, &struc1->cd, mlib->hbond, mlib->hbond_cutoff,
                    2, feat->nbins);
  case 93:
  case 95:
  case 97:
  case 99:
    return itable(f_int1_pt(&struc1->iacc), seq1->nres, ir1, feat->nbins);
  case 94:
  case 96:
  case 98:
  case 100:
    return itable(f_int1_pt(&struc2->iacc), seq2->nres, ir2, feat->nbins);
  case 101:
    tup = property_one_tuple(aln, is1, prop, mlib, ibnd1, ia1, libs);
    return index_inrange(tup->tupclass, feat);
  case 102:
    tup = property_one_tuple(aln, is1, prop, mlib, ibnd1p, ia1p, libs);
    return index_inrange(tup->tupclass, feat);
  case 104:
    if (!tuple_require_natom(mlib, 2, ifi, err)) {
      return 0;
    }
    tup = property_one_tuple(aln, is1, prop, mlib, ibnd1p, ia1p, libs);
    return iangle0(ia1, ia1p, tup->iata[0], struc1, mlib, ifi, feat->nbins);
  case 105:
    if (!tuple_require_natom(mlib, 2, ifi, err)) {
      return 0;
    }
    tup = property_one_tuple(aln, is1, prop, mlib, ibnd1, ia1, libs);
    return iangle0(tup->iata[0], ia1, ia1p, struc1, mlib, ifi, feat->nbins);
  case 106:
    if (!tuple_require_natom(mlib, 2, ifi, err)) {
      return 0;
    }
    tup = property_one_tuple(aln, is1, prop, mlib, ibnd1, ia1, libs);
    tup2 = property_one_tuple(aln, is1, prop, mlib, ibnd1p, ia1p, libs);
    return idihedral0(tup->iata[0], ia1, ia1p, tup2->iata[0], struc1,
                      mlib, ifi, feat->nbins);
  case 107:
    if (!tuple_require_natom(mlib, 3, ifi, err)) {
      return 0;
    }
    tup = property_one_tuple(aln, is1, prop, mlib, ibnd1, ia1, libs);
    tup2 = property_one_tuple(aln, is1, prop, mlib, ibnd1p, ia1p, libs);
    return idihedral0(tup->iata[1], tup2->iata[0], ia1, ia1p, struc1,
                      mlib, ifi, feat->nbins);
  case 108:
    if (!tuple_require_natom(mlib, 3, ifi, err)) {
      return 0;
    }
    tup = property_one_tuple(aln, is1, prop, mlib, ibnd1, ia1, libs);
    tup2 = property_one_tuple(aln, is1, prop, mlib, ibnd1p, ia1p, libs);
    return idihedral0(ia1, ia1p, tup->iata[0], tup2->iata[1], struc1,
                      mlib, ifi, feat->nbins);
  case 109:
    bond = property_one_bond(aln, is1, prop, mlib, MDT_BOND_TYPE_BOND, ibnd1,
                             libs);
    return index_inrange(bond->bndgrp, feat);
  case 110:
    bond = property_one_bond(aln, is1, prop, mlib, MDT_BOND_TYPE_BOND, ibnd1,
                             libs);
    return idist0(bond->iata[0], bond->iata[1], struc1, mlib, ifi,
                  feat->nbins);
  case 111:
    bond = property_one_bond(aln, is1, prop, mlib, MDT_BOND_TYPE_ANGLE, ibnd1,
                             libs);
    return index_inrange(bond->bndgrp, feat);
  case 112:
    bond = property_one_bond(aln, is1, prop, mlib, MDT_BOND_TYPE_ANGLE, ibnd1,
                             libs);
    return iangle0(bond->iata[0], bond->iata[1], bond->iata[2], struc1, mlib,
                   ifi, feat->nbins);
  case 113:
    bond = property_one_bond(aln, is1, prop, mlib, MDT_BOND_TYPE_DIHEDRAL,
                             ibnd1, libs);
    return index_inrange(bond->bndgrp, feat);
  case 114:
    bond = property_one_bond(aln, is1, prop, mlib, MDT_BOND_TYPE_DIHEDRAL,
                             ibnd1, libs);
    return idihedral0(bond->iata[0], bond->iata[1], bond->iata[2],
                      bond->iata[3], struc1, mlib, ifi, feat->nbins);
  case 115:
    irad = property_radius_gyration(aln, is1, prop, mlib, ifi, feat);
    return index_inrange(irad, feat);
  default:
    /* If we don't implement this feature, maybe Modeller does */
    ret = mdt_index(ifi, aln, is1 + 1, ip1 + 1, is2 + 1, ir1 + 1, ir2 + 1,
                    ir1p + 1, ir2p + 1, &mlib->base, ip2 + 1, is3 + 1, ir3 + 1,
                    ir3p + 1, libs, edat, &ierr);
    if (ierr) {
      handle_modeller_error(err);
    }
    return ret;
  }
}
