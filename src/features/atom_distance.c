/** \file atom_distance.c     Atom-atom distance feature.
 *
 *             Part of MDT, Copyright(c) 1989-2008 Andrej Sali
 */

#include "modeller.h"
#include "../mdt_index.h"
#include "../mdt_feature.h"
#include "../geometry.h"

static int getbin(const struct mod_alignment *aln, int protein,
                  int atom1, int atom2, struct mdt_properties *prop, void *data,
                  const struct mod_mdt_libfeature *feat,
                  const struct mdt_library *mlib,
                  const struct mod_libraries *libs, GError **err)
{
  struct mod_structure *s = mod_alignment_structure_get(aln, protein);
  return idist0(atom1, atom2, s, feat);
}

int mdt_feature_atom_distance(struct mdt_library *mlib)
{
  int ifeat;
  ifeat = mdt_feature_atom_pair_add(mlib, "Atom-atom distance",
                                    MOD_MDTC_NONE, FALSE, getbin, NULL);
  mdt_feature_add_needed_file(mlib, ifeat, MOD_MDTF_STRUCTURE);
  return ifeat;
}