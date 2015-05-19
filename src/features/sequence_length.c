/** \file sequence_length.c     Protein sequence length feature.
 *
 *             Part of MDT, Copyright(c) 1989-2015 Andrej Sali
 */

#include "modeller.h"
#include "../mdt_index.h"
#include "../mdt_feature.h"
#include "../mdt_all_features.h"

static int getbin(const struct mod_alignment *aln, int protein,
                  struct mdt_properties *prop,
                  const struct mdt_feature *feat,
                  const struct mdt_library *mlib,
                  const struct mod_libraries *libs, GError **err)
{
  struct mod_sequence *seq = mod_alignment_sequence_get(aln, protein);
  return feat_to_bin(seq->nres, feat);
}

int mdt_feature_sequence_length(struct mdt_library *mlib, int protein,
                                GError **err)
{
  return mdt_feature_protein_add(mlib, "Sequence length", MOD_MDTC_NONE,
                                 protein, getbin, NULL, NULL, err);
}
