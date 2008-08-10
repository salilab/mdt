/** \file sequence_length.c     Protein sequence length feature.
 *
 *             Part of MDT, Copyright(c) 1989-2008 Andrej Sali
 */

#include "modeller.h"
#include "../mdt_index.h"
#include "../mdt_feature.h"

static int getbin(const struct mod_alignment *aln, int protein,
                  struct mdt_properties *prop, void *data,
                  const struct mod_mdt_libfeature *feat,
                  const struct mdt_library *mlib,
                  const struct mod_libraries *libs, GError **err)
{
  struct mod_sequence *seq = mod_alignment_sequence_get(aln, protein);
  return iclsbin(seq->nres, feat);
}

int mdt_feature_sequence_length(struct mdt_library *mlib, int protein,
                                GError **err)
{
  return mdt_feature_protein_add(mlib, "Sequence length", MOD_MDTC_NONE,
                                 protein, getbin, NULL, err);
}
