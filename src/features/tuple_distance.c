/** \file tuple_distance.c     Tuple-tuple non-bonded distance feature.
 *
 *             Part of MDT, Copyright(c) 1989-2020 Andrej Sali
 */

#include "modeller.h"
#include "../mdt_index.h"
#include "../mdt_feature.h"
#include "../mdt_all_features.h"
#include "../mdt_atom_classes.h"
#include "../geometry.h"

static int getbin(const struct mod_alignment *aln, int protein,
                  int atom1, const struct mdt_tuple *tuple1,
                  int atom2, const struct mdt_tuple *tuple2,
                  struct mdt_properties *prop,
                  const struct mdt_feature *feat,
                  const struct mdt_library *mlib,
                  const struct mod_libraries *libs, GError **err)
{
  struct mod_structure *s = mod_alignment_structure_get(aln, protein);
  return idist0(atom1, atom2, s, feat);
}

int mdt_feature_tuple_distance(struct mdt_library *mlib)
{
  int ifeat;
  ifeat = mdt_feature_tuple_pair_add(mlib, "Tuple-tuple non-bonded distance",
                                     MOD_MDTC_NONE, getbin, NULL, NULL);
  mdt_feature_add_needed_file(mlib, ifeat, MOD_MDTF_STRUCTURE);
  return ifeat;
}
