/** \file mdt_index.c  Functions to calculate MDT indices.
 *
 *             Part of MDT, Copyright(c) 1989-2008 Andrej Sali
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

static int irestab(const struct mod_int2_array *ialn, int naln, int iseq,
                   const int *irestyp, int nres, int ip, int delta,
                   gboolean delta_ali, int ndimen, int igaptyp)
{
  if (delta_ali) {
    int ipos = ip + delta;
    if (ipos < 0 || ipos >= naln) {
      return igaptyp;
    } else {
      int ires = mod_int2_get(ialn, ipos, iseq);
      return ires_get(ires, nres, igaptyp, irestyp, ndimen);
    }
  } else {
    int ires = mod_int2_get(ialn, ip, iseq);
    if (ires < 1 || ires > nres) {
      return igaptyp;
    } else {
      return ires_get(ires + delta, nres, igaptyp, irestyp, ndimen);
    }
  }
}

/** Return the bin index itab[ir], or the undefined bin index if anything is
    out of range */
static int itable(const int *itab, int nr, int ir,
                  const struct mod_mdt_libfeature *feat)
{
  if (ir >= 0 && ir < nr && itab[ir] >= 1 && itab[ir] <= feat->nbins) {
    return itab[ir];
  } else {
    return feat->nbins;
  }
}

/** Ensure that a given bin index is in range for the feature; return in
    the undefined bin if not. */
static int index_inrange(int index, const struct mod_mdt_libfeature *feat)
{
  return (index >= 1 && index < feat->nbins) ? index : feat->nbins;
}

/** Convert a raw number to the corresponding feature's MDT bin index */
int iclsbin(float x, const struct mod_mdt_libfeature *feat)
{
  int i;
  const struct mod_mdt_bin *bin = feat->bins;
  for (i = 1; i < feat->nbins; i++, bin++) {
    if (x >= bin->rang1 && x < bin->rang2) {
      return i;
    }
  }
  bin = &feat->bins[0];
  mod_logwarning("iclsbin", "Undefined value; X,x1,x2,n,bin: %f %f %f %d",
                 x, bin->rang1, bin->rang2, feat->nbins);
  return feat->nbins;
}

/** Return the bin index for the raw feature ftab[ir], or the undefined bin
    index if anything is out of range */
int ftable(const float *ftab, int nr, int ir,
           const struct mod_mdt_libfeature *feat)
{
  if (ir >= 0 && ir < nr) {
    return iclsbin(ftab[ir], feat);
  } else {
    return feat->nbins;
  }
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
static int idist0(int ia1, int ia1p, const struct mod_structure *struc,
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

/** Return the bin index for the angle between three specified atoms in the
    same protein. */
static int iangle0(int ia1, int ia2, int ia3,
                   const struct mod_structure *struc,
                   const struct mod_mdt_libfeature *feat)
{
  if (ia1 >= 0 && ia2 >= 0 && ia3 >= 0) {
    float d, *x, *y, *z;
    x = mod_float1_pt(&struc->cd.x);
    y = mod_float1_pt(&struc->cd.y);
    z = mod_float1_pt(&struc->cd.z);
    d = angle1(x[ia1], y[ia1], z[ia1], x[ia2], y[ia2], z[ia2], x[ia3],
               y[ia3], z[ia3]);
    return iclsbin(d, feat);
  } else {
    return feat->nbins;
  }
}

/** Return the bin index for the dihedral angle between four specified atoms
    in the same protein. */
static int idihedral0(int ia1, int ia2, int ia3, int ia4,
                      const struct mod_structure *struc,
                      const struct mod_mdt_libfeature *feat)
{
  if (ia1 >= 0 && ia2 >= 0 && ia3 >= 0 && ia4 >= 0) {
    gboolean outrange;
    float d, *x, *y, *z;
    x = mod_float1_pt(&struc->cd.x);
    y = mod_float1_pt(&struc->cd.y);
    z = mod_float1_pt(&struc->cd.z);
    d = dihedral1(x[ia1], y[ia1], z[ia1], x[ia2], y[ia2], z[ia2], x[ia3],
                  y[ia3], z[ia3], x[ia4], y[ia4], z[ia4], &outrange);
    if (outrange) {
      return feat->nbins;
    } else {
      return iclsbin(d, feat);
    }
  } else {
    return feat->nbins;
  }
}

/** Register our MDT feature types */
void mdt_register_features(struct mod_mdt_library *mlib)
{
  mod_mdt_libfeature_register(mlib, 66, "RESIDUE TYPE AT DELTA I IN A (66)",
                              MOD_MDTC_NONE, MOD_MDTP_A, MOD_MDTS_RESIDUE,
                              FALSE, 0);
  mod_mdt_libfeature_register(mlib, 67, "RESIDUE TYPE AT DELTA I IN B (67)",
                              MOD_MDTC_NONE, MOD_MDTP_B, MOD_MDTS_RESIDUE,
                              FALSE, 0);
  mod_mdt_libfeature_register(mlib, 77, "RESIDUE TYPE AT DELTA J IN A (77)",
                              MOD_MDTC_NONE, MOD_MDTP_A, MOD_MDTS_RESIDUE,
                              FALSE, 0);
  mod_mdt_libfeature_register(mlib, 78, "RESIDUE TYPE AT DELTA J IN B (78)",
                              MOD_MDTC_NONE, MOD_MDTP_B, MOD_MDTS_RESIDUE,
                              FALSE, 0);
  mod_mdt_libfeature_register(mlib, 79, "MODELLER ATOM TYPE OF A (79)",
                              MOD_MDTC_NONE, MOD_MDTP_A, MOD_MDTS_ATOM,
                              FALSE, MOD_MDTF_STRUCTURE, 0);
  mod_mdt_libfeature_register(mlib, 82, "ANY ATOM DISTANCE IN A (82)",
                              MOD_MDTC_NONE, MOD_MDTP_A, MOD_MDTS_ATOM_PAIR,
                              FALSE, MOD_MDTF_STRUCTURE, 0);
  mod_mdt_libfeature_register(mlib, 83, "MODELLER ATOM TYPE AT POS2 OF A (83)",
                              MOD_MDTC_NONE, MOD_MDTP_A, MOD_MDTS_ATOM_PAIR,
                              FALSE, MOD_MDTF_STRUCTURE, 0);
  mod_mdt_libfeature_register(mlib, 84, "H-BOND DONOR IN A (84)",
                              MOD_MDTC_NONE, MOD_MDTP_A, MOD_MDTS_ATOM,
                              FALSE, MOD_MDTF_STRUCTURE, 0);
  mod_mdt_libfeature_register(mlib, 85, "H-BOND ACCEPTOR IN A (85)",
                              MOD_MDTC_NONE, MOD_MDTP_A, MOD_MDTS_ATOM,
                              FALSE, MOD_MDTF_STRUCTURE, 0);
  mod_mdt_libfeature_register(mlib, 86,
                              "HBOND SATISFACTION INDEX OF PROTEIN 1 (86)",
                              MOD_MDTC_NONE, MOD_MDTP_A, MOD_MDTS_PROTEIN,
                              FALSE, MOD_MDTF_STRUCTURE, 0);
  mod_mdt_libfeature_register(mlib, 87, "TOTAL CHARGE AROUND ATOM IN A (87)",
                              MOD_MDTC_NONE, MOD_MDTP_A, MOD_MDTS_ATOM,
                              FALSE, MOD_MDTF_STRUCTURE, 0);
  mod_mdt_libfeature_register(mlib, 101, "ATOM TUPLE TYPE IN A (101)",
                              MOD_MDTC_NONE, MOD_MDTP_A, MOD_MDTS_TUPLE,
                              FALSE, MOD_MDTF_STRUCTURE, 0);
  mod_mdt_libfeature_register(mlib, 102, "ATOM TUPLE TYPE IN A AT POS2 (102)",
                              MOD_MDTC_NONE, MOD_MDTP_A, MOD_MDTS_TUPLE_PAIR,
                              TRUE, MOD_MDTF_STRUCTURE, 0);
  mod_mdt_libfeature_register(mlib, 103,
                              "TUPLE NON-BONDED DISTANCE IN A (103)",
                              MOD_MDTC_NONE, MOD_MDTP_A, MOD_MDTS_TUPLE_PAIR,
                              TRUE, MOD_MDTF_STRUCTURE, 0);
  mod_mdt_libfeature_register(mlib, 104, "NON-BONDED TUPLE ANGLE1 IN A (104)",
                              MOD_MDTC_NONE, MOD_MDTP_A, MOD_MDTS_TUPLE_PAIR,
                              TRUE, MOD_MDTF_STRUCTURE, 0);
  mod_mdt_libfeature_register(mlib, 105, "NON-BONDED TUPLE ANGLE1 IN A (105)",
                              MOD_MDTC_NONE, MOD_MDTP_A, MOD_MDTS_TUPLE_PAIR,
                              TRUE, MOD_MDTF_STRUCTURE, 0);
  mod_mdt_libfeature_register(mlib, 106,
                              "NON-BONDED TUPLE DIHEDRAL1 IN A (106)",
                              MOD_MDTC_NONE, MOD_MDTP_A, MOD_MDTS_TUPLE_PAIR,
                              TRUE, MOD_MDTF_STRUCTURE, 0);
  mod_mdt_libfeature_register(mlib, 107,
                              "NON-BONDED TUPLE DIHEDRAL2 IN A (107)",
                              MOD_MDTC_NONE, MOD_MDTP_A, MOD_MDTS_TUPLE_PAIR,
                              TRUE, MOD_MDTF_STRUCTURE, 0);
  mod_mdt_libfeature_register(mlib, 108,
                              "NON-BONDED TUPLE DIHEDRAL3 IN A (108)",
                              MOD_MDTC_NONE, MOD_MDTP_A, MOD_MDTS_TUPLE_PAIR,
                              TRUE, MOD_MDTF_STRUCTURE, 0);
  mod_mdt_libfeature_register(mlib, 109, "BOND TYPE IN A (109)",
                              MOD_MDTC_NONE, MOD_MDTP_A, MOD_MDTS_BOND,
                              FALSE, MOD_MDTF_STRUCTURE, 0);
  mod_mdt_libfeature_register(mlib, 110, "BOND LENGTH IN A (110)",
                              MOD_MDTC_NONE, MOD_MDTP_A, MOD_MDTS_BOND,
                              FALSE, MOD_MDTF_STRUCTURE, 0);
  mod_mdt_libfeature_register(mlib, 111, "BOND ANGLE TYPE IN A (111)",
                              MOD_MDTC_NONE, MOD_MDTP_A, MOD_MDTS_ANGLE,
                              FALSE, MOD_MDTF_STRUCTURE, 0);
  mod_mdt_libfeature_register(mlib, 112, "BOND ANGLE IN A (112)",
                              MOD_MDTC_NONE, MOD_MDTP_A, MOD_MDTS_ANGLE,
                              FALSE, MOD_MDTF_STRUCTURE, 0);
  mod_mdt_libfeature_register(mlib, 113, "BOND DIHEDRAL ANGLE TYPE IN A (113)",
                              MOD_MDTC_NONE, MOD_MDTP_A, MOD_MDTS_DIHEDRAL,
                              FALSE, MOD_MDTF_STRUCTURE, 0);
  mod_mdt_libfeature_register(mlib, 114, "BOND DIHEDRAL ANGLE IN A (114)",
                              MOD_MDTC_NONE, MOD_MDTP_A, MOD_MDTS_DIHEDRAL,
                              FALSE, MOD_MDTF_STRUCTURE, 0);
}

/** Get the index into the MDT for the given alignment feature */
int my_mdt_index(int ifi, const struct mod_alignment *aln, int is1, int ip1,
                 int is2, int ir1, int ir2, int ir1p, int ir2p, int ia1,
                 int ia1p, const struct mdt_library *mlib, int ip2,
                 int ibnd1, int ibnd1p, int is3, int ir3, int ir3p,
                 const struct mod_libraries *libs,
                 const struct mod_energy_data *edat,
                 struct mdt_properties *prop, GError **err)
{
  int ret, ierr = 0;
  const int *binprop;
  int ibin, ires, iseq, nres;
  float fprop;
  struct mod_structure *struc1, *struc2;
  struct mod_sequence *seq1, *seq2;
  const struct mdt_bond *bond;
  const struct mdt_tuple *tup, *tup2;
  struct mod_mdt_libfeature *feat = &mlib->base.features[ifi - 1];
  struct mdt_feature *mfeat = &g_array_index(mlib->features,
                                             struct mdt_feature, ifi - 1);
  struc1 = mod_alignment_structure_get(aln, is1);
  struc2 = mod_alignment_structure_get(aln, is2);
  seq1 = mod_alignment_sequence_get(aln, is1);
  seq2 = mod_alignment_sequence_get(aln, is2);
  switch (mfeat->type) {
  case MDT_FEATURE_NONE:
    break;
  case MDT_FEATURE_PROTEIN:
    ibin = mfeat->u.protein.getbin(aln, feat->iknown == MOD_MDTP_A ? is1 : is2,
                                   prop, mfeat->data, feat, libs, err);
    if (ibin < 0) {
      return -1;
    } else {
      return index_inrange(ibin, feat);
    }
  case MDT_FEATURE_RESIDUE:
    ires = mfeat->u.residue.delta;
    if (feat->iknown == MOD_MDTP_A) {
      iseq = is1;
      nres = seq1->nres;
      ires += mfeat->u.residue.pos2 ? ir1p : ir1;
    } else {
      iseq = is2;
      nres = seq2->nres;
      ires += mfeat->u.residue.pos2 ? ir2p : ir2;
    }
    if (ires < 0 || ires >= nres) {
      return feat->nbins;
    } else {
      ibin = mfeat->u.residue.getbin(aln, iseq, ires, prop, mfeat->data, feat,
                                     libs, err);
      if (ibin < 0) {
        return -1;
      } else {
        return index_inrange(ibin, feat);
      }
    }
  case MDT_FEATURE_ATOM:
    ibin = mfeat->u.atom.getbin(aln, is1, mfeat->u.atom.pos2 ? ia1p : ia1,
                                prop, mfeat->data, feat, libs, err);
    if (ibin < 0) {
      return -1;
    } else {
      return index_inrange(ibin, feat);
    }
  }
  switch (ifi) {
  case 66:
    return irestab(&aln->ialn, aln->naln, is1, mod_int1_pt(&seq1->irestyp),
                   seq1->nres, ip1, mlib->deltai, mlib->deltai_ali,
                   feat->nbins, libs->igaptyp);
  case 67:
    return irestab(&aln->ialn, aln->naln, is2, mod_int1_pt(&seq2->irestyp),
                   seq2->nres, ip1, mlib->deltai, mlib->deltai_ali,
                   feat->nbins, libs->igaptyp);
  case 77:
    return irestab(&aln->ialn, aln->naln, is1, mod_int1_pt(&seq1->irestyp),
                   seq1->nres, ip1, mlib->deltaj, mlib->deltaj_ali,
                   feat->nbins, libs->igaptyp);
  case 78:
    return irestab(&aln->ialn, aln->naln, is2, mod_int1_pt(&seq2->irestyp),
                   seq2->nres, ip1, mlib->deltaj, mlib->deltaj_ali,
                   feat->nbins, libs->igaptyp);
  case 79:
    binprop = property_iatta(aln, is1, prop, mlib, libs, err);
    if (!binprop) {
      return 0.;
    }
    return itable(binprop, struc1->cd.natm, ia1, feat);
  case 82:
  case 103:
    return idist0(ia1, ia1p, struc1, feat);
  case 83:
    binprop = property_iatta(aln, is1, prop, mlib, libs, err);
    if (!binprop) {
      return 0;
    }
    return itable(binprop, struc1->cd.natm, ia1p, feat);
  case 84:
    binprop = property_hb_iatta(aln, is1, prop, mlib, libs, err);
    if (!binprop) {
      return 0;
    }
    return numb_hda(ia1, binprop, &struc1->cd, mlib->hbond, mlib->hbond_cutoff,
                    0, feat->nbins);
  case 85:
    binprop = property_hb_iatta(aln, is1, prop, mlib, libs, err);
    if (!binprop) {
      return 0;
    }
    return numb_hda(ia1, binprop, &struc1->cd, mlib->hbond, mlib->hbond_cutoff,
                    1, feat->nbins);
  case 86:
    if (!property_hbpot(aln, is1, prop, mlib, libs, &fprop, err)) {
      return 0;
    }
    return iclsbin(fprop, feat);
  case 87:
    binprop = property_hb_iatta(aln, is1, prop, mlib, libs, err);
    if (!binprop) {
      return 0;
    }
    return numb_hda(ia1, binprop, &struc1->cd, mlib->hbond, mlib->hbond_cutoff,
                    2, feat->nbins);
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
    return iangle0(ia1, ia1p, tup->iata[0], struc1, feat);
  case 105:
    if (!tuple_require_natom(mlib, 2, ifi, err)) {
      return 0;
    }
    tup = property_one_tuple(aln, is1, prop, mlib, ibnd1, ia1, libs);
    return iangle0(tup->iata[0], ia1, ia1p, struc1, feat);
  case 106:
    if (!tuple_require_natom(mlib, 2, ifi, err)) {
      return 0;
    }
    tup = property_one_tuple(aln, is1, prop, mlib, ibnd1, ia1, libs);
    tup2 = property_one_tuple(aln, is1, prop, mlib, ibnd1p, ia1p, libs);
    return idihedral0(tup->iata[0], ia1, ia1p, tup2->iata[0], struc1, feat);
  case 107:
    if (!tuple_require_natom(mlib, 3, ifi, err)) {
      return 0;
    }
    tup = property_one_tuple(aln, is1, prop, mlib, ibnd1, ia1, libs);
    tup2 = property_one_tuple(aln, is1, prop, mlib, ibnd1p, ia1p, libs);
    return idihedral0(tup->iata[1], tup2->iata[0], ia1, ia1p, struc1, feat);
  case 108:
    if (!tuple_require_natom(mlib, 3, ifi, err)) {
      return 0;
    }
    tup = property_one_tuple(aln, is1, prop, mlib, ibnd1, ia1, libs);
    tup2 = property_one_tuple(aln, is1, prop, mlib, ibnd1p, ia1p, libs);
    return idihedral0(ia1, ia1p, tup->iata[0], tup2->iata[1], struc1, feat);
  case 109:
    bond = property_one_bond(aln, is1, prop, mlib, MDT_BOND_TYPE_BOND, ibnd1,
                             libs);
    return index_inrange(bond->bndgrp, feat);
  case 110:
    bond = property_one_bond(aln, is1, prop, mlib, MDT_BOND_TYPE_BOND, ibnd1,
                             libs);
    return idist0(bond->iata[0], bond->iata[1], struc1, feat);
  case 111:
    bond = property_one_bond(aln, is1, prop, mlib, MDT_BOND_TYPE_ANGLE, ibnd1,
                             libs);
    return index_inrange(bond->bndgrp, feat);
  case 112:
    bond = property_one_bond(aln, is1, prop, mlib, MDT_BOND_TYPE_ANGLE, ibnd1,
                             libs);
    return iangle0(bond->iata[0], bond->iata[1], bond->iata[2], struc1, feat);
  case 113:
    bond = property_one_bond(aln, is1, prop, mlib, MDT_BOND_TYPE_DIHEDRAL,
                             ibnd1, libs);
    return index_inrange(bond->bndgrp, feat);
  case 114:
    bond = property_one_bond(aln, is1, prop, mlib, MDT_BOND_TYPE_DIHEDRAL,
                             ibnd1, libs);
    return idihedral0(bond->iata[0], bond->iata[1], bond->iata[2],
                      bond->iata[3], struc1, feat);
  default:
    /* If we don't implement this feature, maybe Modeller does */
    ret = mod_mdt_index(ifi, aln, is1 + 1, ip1 + 1, is2 + 1, ir1 + 1, ir2 + 1,
                        ir1p + 1, ir2p + 1, &mlib->base, ip2 + 1, is3 + 1,
                        ir3 + 1, ir3p + 1, libs, edat, &ierr);
    if (ierr) {
      handle_modeller_error(err);
    }
    return ret;
  }
}
