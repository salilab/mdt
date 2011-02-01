/** \file average_residue_accessibility.c Average residue accessibility feature.
 *
 *             Part of MDT, Copyright(c) 1989-2011 Andrej Sali
 */

#include "modeller.h"
#include "../mdt_index.h"
#include "../mdt_feature.h"
#include "../mdt_all_features.h"

static int getbin(const struct mod_alignment *aln, int protein, int residue1,
                  int residue2, struct mdt_properties *prop, void *data,
                  const struct mod_mdt_libfeature *feat,
                  const struct mdt_library *mlib,
                  const struct mod_libraries *libs, GError **err)
{
  struct mod_structure *s = mod_alignment_structure_get(aln, protein);
  float f1 = mod_float1_get(&s->acc, residue1);
  float f2 = mod_float1_get(&s->acc, residue2);
  return feat_to_bin(0.5 * (f1 + f2), feat);
}

int mdt_feature_average_residue_accessibility(struct mdt_library *mlib,
                                              int protein, GError **err)
{
  int ifeat;
  ifeat = mdt_feature_residue_pair_add(
              mlib, "Average accessibility of a residue pair", MOD_MDTC_NONE,
              protein, FALSE, getbin, NULL, NULL, err);
  if (ifeat >= 0) {
    mdt_feature_add_needed_file(mlib, ifeat, MOD_MDTF_STRUCTURE);
    mdt_feature_add_needed_file(mlib, ifeat, MOD_MDTF_PSA);
  }
  return ifeat;
}
