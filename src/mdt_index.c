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
  int ibin, ires, iseq, nres, iatom, ibnd, ialnpos, ires1, ires2;
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
  case MDT_FEATURE_PROTEIN_PAIR:
    ibin = mfeat->u.protein_pair.getbin(aln, is1,
                                        feat->iknown == MOD_MDTP_AB ? is2 : is3,
                                        prop, mfeat->data, feat, mlib, libs,
                                        err);
    if (ibin < 0) {
      return -1;
    } else {
      return index_inrange(ibin, feat);
    }
  case MDT_FEATURE_RESIDUE:
    if (feat->iknown == MOD_MDTP_A) {
      iseq = is1;
      nres = seq1->nres;
      ires = mfeat->u.residue.pos2 ? ir1p : ir1;
    } else {
      iseq = is2;
      nres = seq2->nres;
      ires = mfeat->u.residue.pos2 ? ir2p : ir2;
    }
    ires += mfeat->u.residue.delta;
    if (ires < 0 || ires >= nres) {
      return index_inrange(mfeat->u.residue.bin_seq_outrange, feat);
    }
    if (mfeat->u.residue.align_delta != 0) {
      /* Convert to alignment position, and apply align_delta */
      ialnpos = mod_int2_get(&aln->invaln, ires, iseq) - 1;
      ialnpos += mfeat->u.residue.align_delta;
      if (ialnpos < 0 || ialnpos >= aln->naln) {
        return index_inrange(mfeat->u.residue.bin_seq_outrange, feat);
      }
      /* Convert back to residue position */
      ires = mod_int2_get(&aln->ialn, ialnpos, iseq) - 1;
      if (ires < 0 || ires >= nres) {
        return index_inrange(mfeat->u.residue.bin_seq_outrange, feat);
      }
    }
    ibin = mfeat->u.residue.getbin(aln, iseq, ires, prop,
                                   mfeat->data, feat, mlib, libs, err);
    if (ibin < 0) {
      return -1;
    } else {
      return index_inrange(ibin, feat);
    }
  case MDT_FEATURE_ALIGNED_RESIDUE:
    if (feat->iknown == MOD_MDTP_AB) {
      iseq = is2;
    } else {
      iseq = is3;
    }
    ibin = mfeat->u.aligned_residue.getbin(aln, is1, iseq, ip1, prop,
                                           mfeat->data, feat, mlib, libs, err);
    if (ibin < 0) {
      return -1;
    } else {
      return index_inrange(ibin, feat);
    }
  case MDT_FEATURE_RESIDUE_PAIR:
    if (feat->iknown == MOD_MDTP_A) {
      iseq = is1;
      ires1 = ir1;
      ires2 = ir1p;
    } else {
      iseq = is2;
      ires1 = ir2;
      ires2 = ir2p;
    }
    ibin = mfeat->u.residue_pair.getbin(aln, iseq, ires1, ires2, prop,
                                        mfeat->data, feat, mlib, libs, err);
    if (ibin < 0) {
      return -1;
    } else {
      return index_inrange(ibin, feat);
    }
  case MDT_FEATURE_ALIGNED_RESIDUE_PAIR:
    if (feat->iknown == MOD_MDTP_AB) {
      iseq = is2;
    } else {
      iseq = is3;
    }
    ibin = mfeat->u.aligned_residue_pair.getbin(aln, is1, iseq, ip1, ip2, prop,
                                                mfeat->data, feat, mlib, libs,
                                                err);
    if (ibin < 0) {
      return -1;
    } else {
      return index_inrange(ibin, feat);
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
  /* If we don't implement this feature, maybe Modeller does */
  ret = mod_mdt_index(ifi, aln, is1 + 1, ip1 + 1, is2 + 1, ir1 + 1, ir2 + 1,
                      ir1p + 1, ir2p + 1, &mlib->base, ip2 + 1, is3 + 1,
                      ir3 + 1, ir3p + 1, libs, edat, &ierr);
  if (ierr) {
    handle_modeller_error(err);
  }
  return ret;
}
