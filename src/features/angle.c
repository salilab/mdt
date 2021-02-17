/** \file angle.c  Angle feature.
 *
 *             Part of MDT, Copyright(c) 1989-2021 Andrej Sali
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
  return iangle0(bond->iata[0], bond->iata[1], bond->iata[2], s, feat);
}

int mdt_feature_angle(struct mdt_library *mlib)
{
  int ifeat;
  ifeat = mdt_feature_angle_add(mlib, "Angle", MOD_MDTC_NONE, getbin, NULL,
                                NULL);
  mdt_feature_add_needed_file(mlib, ifeat, MOD_MDTF_STRUCTURE);
  return ifeat;
}
