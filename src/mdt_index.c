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
#include "mdt_triplets.h"
#include "mdt_geometry.h"

static int ires_get(int ires, int nres, int igaptyp, const int *irestyp,
                    int ndimen)
{
  if (ires < 1 || ires > nres) {
    return igaptyp;
  } else if (irestyp[ires - 1] >= 21 || irestyp[ires - 1] <= 0) {
    modlogwarning("irestab", "ires, irestp: %d %d", ires, irestyp[ires - 1]);
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

static int itable(const int *itab, int nr, int ir, int ndim)
{
  if (ir >= 0 && ir < nr && itab[ir] >= 1 && itab[ir] <= ndim) {
    return itab[ir];
  } else {
    return ndim;
  }
}

/** Convert a raw number to the corresponding MDT bin index */
int iclsbin(float x, const struct mdt_library *mlib, int ifi, int nrang)
{
  int i;
  const struct mdt_libfeature *feat = &mlib->base.features[ifi - 1];
  const struct mdt_bin *bin = feat->bins;
  for (i = 0; i < nrang; i++, bin++) {
    if (x >= bin->rang1 && x <= bin->rang2) {
      return i + 1;
    }
  }
  bin = &feat->bins[0];
  modlogwarning("iclsbin", "Undefined value; X,x1,x2,n,bin: %f %f %f %d",
                x, bin->rang1, bin->rang2, nrang + 1);
  return nrang + 1;
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
    return nrang + 1;
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
    return nrang + 1;
  }
}

/** Return the bin index for the dihedral angle between four specified atoms
    in the same protein. */
static int idihedral0(int ia1, int ia2, int ia3, int ia4,
                      const struct structure *struc,
                      const struct mdt_library *mlib, int ifi, int nrang)
{
  if (ia1 >= 0 && ia2 >= 0 && ia3 >= 0 && ia4 >= 0) {
    float d, *x, *y, *z;
    x = f_float1_pt(&struc->cd.x);
    y = f_float1_pt(&struc->cd.y);
    z = f_float1_pt(&struc->cd.z);
    d = dihedral1(x[ia1], y[ia1], z[ia1], x[ia2], y[ia2], z[ia2], x[ia3],
                  y[ia3], z[ia3], x[ia4], y[ia4], z[ia4]);
    return iclsbin(d, mlib, ifi, nrang);
  } else {
    return nrang + 1;
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
  int iresol;
  float fprop;
  struct structure *struc1, *struc2;
  struct sequence *seq1, *seq2;
  const struct mdt_bond *bond;
  const struct mdt_triplet *trp, *trp2;
  struct mdt_libfeature *feat = &mlib->base.features[ifi - 1];
  struc1 = alignment_structure_get(aln, is1);
  struc2 = alignment_structure_get(aln, is2);
  seq1 = alignment_sequence_get(aln, is1);
  seq2 = alignment_sequence_get(aln, is2);
  switch (ifi) {
  case 35:
    iresol = property_iresol(aln, is1, prop, mlib, ifi, feat);
    return itable(&iresol, 1, 0, feat->nbins);
  case 38:
    iresol = property_iresol(aln, is2, prop, mlib, ifi, feat);
    return itable(&iresol, 1, 0, feat->nbins);
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
    return iclsbin(fprop, mlib, ifi, feat->nbins - 1);
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
    trp = property_one_triplet(aln, is1, prop, mlib, ibnd1, ia1, libs);
    return CLAMP(trp->trpclass, 1, feat->nbins);
  case 102:
    trp = property_one_triplet(aln, is1, prop, mlib, ibnd1p, ia1p, libs);
    return CLAMP(trp->trpclass, 1, feat->nbins);
  case 104:
    trp = property_one_triplet(aln, is1, prop, mlib, ibnd1p, ia1p, libs);
    return iangle0(ia1, ia1p, trp->iata[0], struc1, mlib, ifi, feat->nbins);
  case 105:
    trp = property_one_triplet(aln, is1, prop, mlib, ibnd1, ia1, libs);
    return iangle0(trp->iata[0], ia1, ia1p, struc1, mlib, ifi, feat->nbins);
  case 106:
    trp = property_one_triplet(aln, is1, prop, mlib, ibnd1, ia1, libs);
    trp2 = property_one_triplet(aln, is1, prop, mlib, ibnd1p, ia1p, libs);
    return idihedral0(trp->iata[0], ia1, ia1p, trp2->iata[0], struc1,
                      mlib, ifi, feat->nbins);
  case 107:
    trp = property_one_triplet(aln, is1, prop, mlib, ibnd1, ia1, libs);
    trp2 = property_one_triplet(aln, is1, prop, mlib, ibnd1p, ia1p, libs);
    return idihedral0(trp->iata[1], trp2->iata[0], ia1, ia1p, struc1,
                      mlib, ifi, feat->nbins);
  case 108:
    trp = property_one_triplet(aln, is1, prop, mlib, ibnd1, ia1, libs);
    trp2 = property_one_triplet(aln, is1, prop, mlib, ibnd1p, ia1p, libs);
    return idihedral0(ia1, ia1p, trp->iata[0], trp2->iata[1], struc1,
                      mlib, ifi, feat->nbins);
  case 109:
    bond = property_one_bond(aln, is1, prop, mlib, MDT_BOND_TYPE_BOND, ibnd1,
                             libs);
    return CLAMP(bond->bndgrp, 1, feat->nbins);
  case 110:
    bond = property_one_bond(aln, is1, prop, mlib, MDT_BOND_TYPE_BOND, ibnd1,
                             libs);
    return idist0(bond->iata[0], bond->iata[1], struc1, mlib, ifi,
                  feat->nbins);
  case 111:
    bond = property_one_bond(aln, is1, prop, mlib, MDT_BOND_TYPE_ANGLE, ibnd1,
                             libs);
    return CLAMP(bond->bndgrp, 1, feat->nbins);
  case 112:
    bond = property_one_bond(aln, is1, prop, mlib, MDT_BOND_TYPE_ANGLE, ibnd1,
                             libs);
    return iangle0(bond->iata[0], bond->iata[1], bond->iata[2], struc1, mlib,
                   ifi, feat->nbins);
  case 113:
    bond = property_one_bond(aln, is1, prop, mlib, MDT_BOND_TYPE_DIHEDRAL,
                             ibnd1, libs);
    return CLAMP(bond->bndgrp, 1, feat->nbins);
  case 114:
    bond = property_one_bond(aln, is1, prop, mlib, MDT_BOND_TYPE_DIHEDRAL,
                             ibnd1, libs);
    return idihedral0(bond->iata[0], bond->iata[1], bond->iata[2],
                      bond->iata[3], struc1, mlib, ifi, feat->nbins);
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
