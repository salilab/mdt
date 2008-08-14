/** \file residue_distance_difference.c
 *  \brief Residue-residue distance difference feature.
 *
 *             Part of MDT, Copyright(c) 1989-2008 Andrej Sali
 */

#include "modeller.h"
#include "../mdt_index.h"
#include "../mdt_feature.h"
#include "../mdt_all_features.h"
#include "../geometry.h"

/** Get the distance between two alignment positions in the same protein */
static float get_distance(const struct mod_alignment *aln, int protein,
                          int alnpos1, int alnpos2)
{
  struct mod_structure *s = mod_alignment_structure_get(aln, protein);
  int residue1 = mod_int2_get(&aln->ialn, alnpos1, protein) - 1;
  int residue2 = mod_int2_get(&aln->ialn, alnpos2, protein) - 1;
  int atom1 = mod_int1_get(&s->idsta1, residue1) - 1;
  int atom2 = mod_int1_get(&s->idsta1, residue2) - 1;
  float *x = mod_float1_pt(&s->cd.x);
  float *y = mod_float1_pt(&s->cd.y);
  float *z = mod_float1_pt(&s->cd.z);
  return dist1(x[atom1], y[atom1], z[atom1], x[atom2], y[atom2], z[atom2]);
}

static int getbin(const struct mod_alignment *aln, int protein1, int protein2,
                  int alnpos1, int alnpos2, struct mdt_properties *prop,
                  void *data, const struct mod_mdt_libfeature *feat,
                  const struct mdt_library *mlib,
                  const struct mod_libraries *libs, GError **err)
{
  float d1 = get_distance(aln, protein1, alnpos1, alnpos2);
  float d2 = get_distance(aln, protein2, alnpos1, alnpos2);
  return feat_to_bin(d2 - d1, feat);
}

int mdt_feature_residue_distance_difference(struct mdt_library *mlib,
                                            int protein1, int protein2,
                                            GError **err)
{
  int ifeat;
  ifeat = mdt_feature_aligned_residue_pair_add(
              mlib, "Residue-residue distance difference", MOD_MDTC_ATIND,
              protein1, protein2, TRUE, getbin, NULL, err);
  if (ifeat >= 0) {
    mdt_feature_add_needed_file(mlib, ifeat, MOD_MDTF_STRUCTURE);
  }
  return ifeat;
}
