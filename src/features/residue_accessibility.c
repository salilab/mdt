/** \file residue_accessibility.c  Residue accessibility feature.
 *
 *             Part of MDT, Copyright(c) 1989-2025 Andrej Sali
 */

#include "modeller.h"
#include "../mdt_index.h"
#include "../mdt_feature.h"
#include "../mdt_all_features.h"

static int getbin(const struct mod_alignment *aln, int protein, int residue,
                  struct mdt_properties *prop,
                  const struct mdt_feature *feat,
                  const struct mdt_library *mlib,
                  const struct mod_libraries *libs, GError **err)
{
  struct mod_structure *s = mod_alignment_structure_get(aln, protein);
  float f = mod_float1_get(&s->acc, residue);
  return feat_to_bin(f, feat);
}

int mdt_feature_residue_accessibility(struct mdt_library *mlib, int protein,
                                      int delta, int align_delta, gboolean pos2,
                                      GError **err)
{
  int ifeat;
  ifeat = mdt_feature_residue_add(mlib, "Residue accessibility", MOD_MDTC_NONE,
                                  protein, delta, align_delta, pos2, -1, getbin,
                                  NULL, NULL, err);
  if (ifeat >= 0) {
    mdt_feature_add_needed_file(mlib, ifeat, MOD_MDTF_STRUCTURE);
    mdt_feature_add_needed_file(mlib, ifeat, MOD_MDTF_PSA);
  }
  return ifeat;
}
