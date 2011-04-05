/** \file dihedral.c  Dihedral feature.
 *
 *             Part of MDT, Copyright(c) 1989-2011 Andrej Sali
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
  return idihedral0(bond->iata[0], bond->iata[1], bond->iata[2], bond->iata[3],
                    s, feat);
}

int mdt_feature_dihedral(struct mdt_library *mlib)
{
  int ifeat;
  ifeat = mdt_feature_dihedral_add(mlib, "Dihedral", MOD_MDTC_NONE, getbin,
                                   NULL, NULL);
  mdt_feature_add_needed_file(mlib, ifeat, MOD_MDTF_STRUCTURE);
  mdt_feature_periodic_set(mlib, ifeat, TRUE);
  return ifeat;
}
