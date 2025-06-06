/** \file bond_length.c  Bond length feature.
 *
 *             Part of MDT, Copyright(c) 1989-2025 Andrej Sali
 */

#include "modeller.h"
#include "../mdt_index.h"
#include "../mdt_feature.h"
#include "../mdt_all_features.h"
#include "../geometry.h"

static int getbin(const struct mod_alignment *aln, int protein,
                  const struct mdt_bond *bond, struct mdt_properties *prop,
                  const struct mdt_feature *feat,
                  const struct mdt_library *mlib,
                  const struct mod_libraries *libs, GError **err)
{
  struct mod_structure *s = mod_alignment_structure_get(aln, protein);
  return idist0(bond->iata[0], bond->iata[1], s, feat);
}

int mdt_feature_bond_length(struct mdt_library *mlib)
{
  int ifeat;
  ifeat = mdt_feature_bond_add(mlib, "Bond length", MOD_MDTC_NONE, getbin,
                               NULL, NULL);
  if (ifeat >= 0) {
    mdt_feature_add_needed_file(mlib, ifeat, MOD_MDTF_STRUCTURE);
  }
  return ifeat;
}
