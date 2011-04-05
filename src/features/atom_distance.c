/** \file atom_distance.c     Atom-atom distance feature.
 *
 *             Part of MDT, Copyright(c) 1989-2011 Andrej Sali
 */

#include "modeller.h"
#include "../mdt_index.h"
#include "../mdt_feature.h"
#include "../mdt_all_features.h"
#include "../geometry.h"

static int getbin(const struct mod_alignment *aln, int protein,
                  int atom1, int atom2, struct mdt_properties *prop,
                  const struct mdt_feature *feat,
                  const struct mdt_library *mlib,
                  const struct mod_libraries *libs, GError **err)
{
  struct mod_structure *s = mod_alignment_structure_get(aln, protein);
  return idist0(atom1, atom2, s, feat);
}

int mdt_feature_atom_distance(struct mdt_library *mlib)
{
  return mdt_feature_atom_pair_add(mlib, "Atom-atom distance",
                                   MOD_MDTC_NONE, FALSE, getbin, NULL, NULL);
}
