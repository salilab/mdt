/** \file residue_index_difference.c    Residue index difference feature.
 *
 *             Part of MDT, Copyright(c) 1989-2012 Andrej Sali
 */

#include "modeller.h"
#include "../mdt_index.h"
#include "../mdt_feature.h"
#include "../mdt_all_features.h"
#include <stdlib.h>

static int getbin(const struct mod_alignment *aln, int protein, int residue1,
                  int residue2, struct mdt_properties *prop,
                  const struct mdt_feature *feat,
                  const struct mdt_library *mlib,
                  const struct mod_libraries *libs, GError **err)
{
  return feat_to_bin(residue2 - residue1, feat);
}

static int absgetbin(const struct mod_alignment *aln, int protein, int residue1,
                     int residue2, struct mdt_properties *prop,
                     const struct mdt_feature *feat,
                     const struct mdt_library *mlib,
                     const struct mod_libraries *libs, GError **err)
{
  return feat_to_bin(abs(residue2 - residue1), feat);
}

int mdt_feature_residue_index_difference(struct mdt_library *mlib,
                                         int protein, gboolean absolute,
                                         GError **err)
{
  if (absolute) {
    return mdt_feature_residue_pair_add(mlib,
                                        "Absolute Residue index difference",
                                        MOD_MDTC_NONE, protein, FALSE,
                                        absgetbin, NULL, NULL, err);
  } else {
    return mdt_feature_residue_pair_add(mlib, "Residue index difference",
                                        MOD_MDTC_NONE, protein, TRUE, getbin,
                                        NULL, NULL, err);
  }
}
