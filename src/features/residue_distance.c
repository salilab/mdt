/** \file residue_distance.c     Residue-residue distance feature.
 *
 *             Part of MDT, Copyright(c) 1989-2018 Andrej Sali
 */

#include "modeller.h"
#include "../mdt_index.h"
#include "../mdt_feature.h"
#include "../mdt_all_features.h"
#include "../mdt_property.h"
#include "../geometry.h"

static int getbin(const struct mod_alignment *aln, int protein, int residue1,
                  int residue2, struct mdt_properties *prop,
                  const struct mdt_feature *feat,
                  const struct mdt_library *mlib,
                  const struct mod_libraries *libs, GError **err)
{
  const int *dstind1, *dstind2;
  struct mod_structure *s = mod_alignment_structure_get(aln, protein);

  property_distance_atom_indices(aln, protein, prop, mlib, &dstind1, &dstind2);
  return idist0(dstind1[residue1], dstind2[residue2], s, feat);
}

int mdt_feature_residue_distance(struct mdt_library *mlib, int protein,
                                 GError **err)
{
  int ifeat;
  ifeat = mdt_feature_residue_pair_add(mlib, "Residue-residue distance",
                                       MOD_MDTC_NONE, protein, TRUE, getbin,
                                       NULL, NULL, err);
  if (ifeat >= 0) {
    mdt_feature_add_needed_file(mlib, ifeat, MOD_MDTF_STRUCTURE);
  }
  return ifeat;
}
