/** \file mdt_index.c  Functions to calculate MDT indices.
 *
 *             Part of MDT, Copyright(c) 1989-2013 Andrej Sali
 */

#include <string.h>
#include <math.h>
#include "modeller.h"
#include "util.h"
#include "geometry.h"
#include "mdt_feature.h"
#include "mdt_index.h"
#include "mdt_property.h"
#include "mdt_hydrogen_bonds.h"
#include "mdt_stereo.h"
#include "mdt_tuples.h"

/** Ensure that a given bin index is in range for the feature; return in
    the undefined bin if not. If the index is -1, indicating an error,
    preserve it. */
static int index_inrange(int index, const struct mod_mdt_libfeature *feat)
{
  return (index >= 1 && index < feat->nbins) ? index : feat->nbins;
}

/** Call index_inrange, unless index indicates an error (-1), in which case
    preserve it. */
static int index_inrange_err(int index, const struct mod_mdt_libfeature *feat)
{
  return index == -1 ? index : index_inrange(index, feat);
}

/** Convert a raw number to the corresponding feature's MDT bin index */
int feat_to_bin(float x, const struct mdt_feature *feat)
{
  const struct mod_mdt_libfeature *base = feat->base;
  const struct mod_mdt_bin *bin = base->bins;
  if (feat->uniform_bins) {
    int index = (int)((x - bin[0].rang1) * feat->inverse_bin_width);
    if (index >= 0 && index < base->nbins - 1) {
      return index + 1;
    }
  } else {
    int i;
    for (i = 1; i < base->nbins; i++, bin++) {
      if (x >= bin->rang1 && x < bin->rang2) {
        return i;
      }
    }
    bin = &base->bins[0];
  }
  mod_logwarning("feat_to_bin", "Undefined value; X,x1,x2,n,bin: %f %f %f %d",
                 x, bin->rang1, bin->rang2,
                 mdt_feature_undefined_bin_get(feat));
  return mdt_feature_undefined_bin_get(feat);
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
  int ibin, ires, iseq, nres, iatom, ibnd, ialnpos, ires1, ires2, ibin1, ibin2;
  struct mod_sequence *seq1, *seq2, *seq3;
  const struct mdt_bond *bond;
  const struct mdt_tuple *tup, *tup2;
  struct mod_mdt_libfeature *feat = &mlib->base.features[ifi - 1];
  const struct mdt_feature *mfeat = &g_array_index(mlib->features,
                                                   const struct mdt_feature,
                                                   ifi - 1);
  switch (mfeat->type) {
  default:
    g_assert_not_reached();
    return -1;
  case MDT_FEATURE_GROUP:
    ibin1 = my_mdt_index(mfeat->u.group.ifeat1, aln, is1, ip1, is2, ir1, ir2,
                         ir1p, ir2p, ia1, ia1p, mlib, ip2, ibnd1, ibnd1p, is3,
                         ir3, ir3p, libs, edat, prop, err);
    // todo: handle error
    ibin2 = my_mdt_index(mfeat->u.group.ifeat2, aln, is1, ip1, is2, ir1, ir2,
                         ir1p, ir2p, ia1, ia1p, mlib, ip2, ibnd1, ibnd1p, is3,
                         ir3, ir3p, libs, edat, prop, err);
    // todo: handle error
    ibin = mfeat->u.group.getbin(ibin1, ibin2, prop, mfeat, mlib, libs, err);
    return index_inrange_err(ibin, feat);
  case MDT_FEATURE_PROTEIN:
    ibin = mfeat->u.protein.getbin(aln, feat->iknown == MOD_MDTP_A ? is1 :
                                        feat->iknown == MOD_MDTP_B ? is2 : is3,
                                   prop, mfeat, mlib, libs, err);
    return index_inrange_err(ibin, feat);
  case MDT_FEATURE_PROTEIN_PAIR:
    ibin = mfeat->u.protein_pair.getbin(aln, is1,
                                        feat->iknown == MOD_MDTP_AB ? is2 : is3,
                                        prop, mfeat, mlib, libs, err);
    return index_inrange_err(ibin, feat);
  case MDT_FEATURE_RESIDUE:
    seq1 = mod_alignment_sequence_get(aln, is1);
    seq2 = mod_alignment_sequence_get(aln, is2);
    seq3 = mod_alignment_sequence_get(aln, is3);
    switch(feat->iknown) {
    case MOD_MDTP_A:
      iseq = is1;
      nres = seq1->nres;
      ires = mfeat->u.residue.pos2 ? ir1p : ir1;
      break;
    case MOD_MDTP_B:
      iseq = is2;
      nres = seq2->nres;
      ires = mfeat->u.residue.pos2 ? ir2p : ir2;
      break;
    default:
      iseq = is3;
      nres = seq3->nres;
      ires = mfeat->u.residue.pos2 ? ir3p : ir3;
      break;
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
                                   mfeat, mlib, libs, err);
    return index_inrange_err(ibin, feat);
  case MDT_FEATURE_ALIGNED_RESIDUE:
    if (feat->iknown == MOD_MDTP_AB) {
      iseq = is2;
    } else {
      iseq = is3;
    }
    ibin = mfeat->u.aligned_residue.getbin(aln, is1, iseq, ip1, prop,
                                           mfeat, mlib, libs, err);
    return index_inrange_err(ibin, feat);
  case MDT_FEATURE_RESIDUE_PAIR:
    switch(feat->iknown) {
    case MOD_MDTP_A:
      iseq = is1;
      ires1 = ir1;
      ires2 = ir1p;
      break;
    case MOD_MDTP_B:
      iseq = is2;
      ires1 = ir2;
      ires2 = ir2p;
      break;
    default:
      iseq = is3;
      ires1 = ir3;
      ires2 = ir3p;
      break;
    }
    ibin = mfeat->u.residue_pair.getbin(aln, iseq, ires1, ires2, prop,
                                        mfeat, mlib, libs, err);
    return index_inrange_err(ibin, feat);
  case MDT_FEATURE_ALIGNED_RESIDUE_PAIR:
    if (feat->iknown == MOD_MDTP_AB) {
      iseq = is2;
    } else {
      iseq = is3;
    }
    ibin = mfeat->u.aligned_residue_pair.getbin(aln, is1, iseq, ip1, ip2, prop,
                                                mfeat, mlib, libs, err);
    return index_inrange_err(ibin, feat);
  case MDT_FEATURE_ATOM:
    ibin = mfeat->u.atom.getbin(aln, is1, mfeat->u.atom.pos2 ? ia1p : ia1,
                                prop, mfeat, mlib, libs, err);
    return index_inrange_err(ibin, feat);
  case MDT_FEATURE_ATOM_PAIR:
    ibin = mfeat->u.atom_pair.getbin(aln, is1, ia1, ia1p, prop, mfeat,
                                     mlib, libs, err);
    return index_inrange_err(ibin, feat);
  case MDT_FEATURE_TUPLE:
    if (mfeat->u.tuple.pos2) {
      iatom = ia1p;
      ibnd = ibnd1p;
    } else {
      iatom = ia1;
      ibnd = ibnd1;
    }
    tup = property_one_tuple(aln, is1, prop, mlib, ibnd, iatom, libs);
    ibin = mfeat->u.tuple.getbin(aln, is1, iatom, tup, prop, mfeat, mlib,
                                 libs, err);
    return index_inrange_err(ibin, feat);
  case MDT_FEATURE_TUPLE_PAIR:
    tup = property_one_tuple(aln, is1, prop, mlib, ibnd1, ia1, libs);
    tup2 = property_one_tuple(aln, is1, prop, mlib, ibnd1p, ia1p, libs);
    ibin = mfeat->u.tuple_pair.getbin(aln, is1, ia1, tup, ia1p, tup2, prop,
                                      mfeat, mlib, libs, err);
    return index_inrange_err(ibin, feat);
  case MDT_FEATURE_BOND:
    bond = property_one_bond(aln, is1, prop, mlib, mfeat->u.bond.type,
                             ibnd1, libs);
    ibin = mfeat->u.bond.getbin(aln, is1, bond, prop, mfeat, mlib, libs, err);
    return index_inrange_err(ibin, feat);
  }
}
