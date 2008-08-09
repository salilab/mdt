/** \file mdt_index.c  Functions to calculate MDT indices.
 *
 *             Part of MDT, Copyright(c) 1989-2008 Andrej Sali
 */

#include <string.h>
#include <math.h>
#include "modeller.h"
#include "util.h"
#include "geometry.h"
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
  int ibin, ires, iseq, nres, iatom, ibnd;
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
                                   prop, mfeat->data, feat, mlib, libs, err);
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
                                     mlib, libs, err);
      if (ibin < 0) {
        return -1;
      } else {
        return index_inrange(ibin, feat);
      }
    }
  case MDT_FEATURE_ATOM:
    ibin = mfeat->u.atom.getbin(aln, is1, mfeat->u.atom.pos2 ? ia1p : ia1,
                                prop, mfeat->data, feat, mlib, libs, err);
    if (ibin < 0) {
      return -1;
    } else {
      return index_inrange(ibin, feat);
    }
  case MDT_FEATURE_ATOM_PAIR:
    ibin = mfeat->u.atom_pair.getbin(aln, is1, ia1, ia1p, prop, mfeat->data,
                                     feat, mlib, libs, err);
    if (ibin < 0) {
      return -1;
    } else {
      return index_inrange(ibin, feat);
    }
  case MDT_FEATURE_TUPLE:
    if (mfeat->u.tuple.pos2) {
      iatom = ia1p;
      ibnd = ibnd1p;
    } else {
      iatom = ia1;
      ibnd = ibnd1;
    }
    tup = property_one_tuple(aln, is1, prop, mlib, ibnd, iatom, libs);
    ibin = mfeat->u.tuple.getbin(aln, is1, iatom, tup, prop, mfeat->data,
                                 feat, mlib, libs, err);
    if (ibin < 0) {
      return -1;
    } else {
      return index_inrange(ibin, feat);
    }
  case MDT_FEATURE_TUPLE_PAIR:
    tup = property_one_tuple(aln, is1, prop, mlib, ibnd1, ia1, libs);
    tup2 = property_one_tuple(aln, is1, prop, mlib, ibnd1p, ia1p, libs);
    ibin = mfeat->u.tuple_pair.getbin(aln, is1, ia1, tup, ia1p, tup2, prop,
                                      mfeat->data, feat, mlib, libs, err);
    if (ibin < 0) {
      return -1;
    } else {
      return index_inrange(ibin, feat);
    }
  case MDT_FEATURE_BOND:
    bond = property_one_bond(aln, is1, prop, mlib, mfeat->u.bond.type,
                             ibnd1, libs);
    ibin = mfeat->u.bond.getbin(aln, is1, bond, prop, mfeat->data, feat, mlib,
                                libs, err);
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
