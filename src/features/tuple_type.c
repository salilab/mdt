/** \file tuple_type.c  Tuple type feature.
 *
 *             Part of MDT, Copyright(c) 1989-2020 Andrej Sali
 */

#include "modeller.h"
#include "../mdt_index.h"
#include "../mdt_feature.h"
#include "../mdt_all_features.h"
#include "../mdt_atom_classes.h"
#include "../mdt_tuples.h"

static int getbin(const struct mod_alignment *aln, int protein, int atom,
                  const struct mdt_tuple *tuple, struct mdt_properties *prop,
                  const struct mdt_feature *feat,
                  const struct mdt_library *mlib,
                  const struct mod_libraries *libs, GError **err)
{
  return tuple->tupclass;
}

int mdt_feature_tuple_type(struct mdt_library *mlib, gboolean pos2)
{
  int ifeat;
  struct mod_mdt_libfeature *feat;
  ifeat = mdt_feature_tuple_add(mlib, "Tuple type", MOD_MDTC_NONE,
                                pos2, getbin, NULL, NULL);
  mdt_feature_add_needed_file(mlib, ifeat, MOD_MDTF_STRUCTURE);
  /* Set number of bins and their symbols */
  feat = &mlib->base.features[ifeat - 1];
  update_mdt_feat_atclass(feat, mlib->tupclass);
  return ifeat;
}
