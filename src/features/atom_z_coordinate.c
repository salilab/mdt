/** \file atom_z_coordinate.c  Z-coordinate of an atom.
 *
 *             Part of MDT, Copyright(c) 1989-2011 Andrej Sali
 */

#include "modeller.h"
#include "../mdt_index.h"
#include "../mdt_feature.h"
#include "../mdt_all_features.h"

static int getbin(const struct mod_alignment *aln, int protein, int atom,
                  struct mdt_properties *prop, void *data,
                  const struct mod_mdt_libfeature *feat,
                  const struct mdt_library *mlib,
                  const struct mod_libraries *libs, GError **err)
{
  struct mod_structure *s = mod_alignment_structure_get(aln, protein);
  float f = mod_float1_get(&s->cd.z, atom);
  return feat_to_bin(f, feat);
}

int mdt_feature_z_coordinate(struct mdt_library *mlib, gboolean pos2)
{
  int ifeat;
  ifeat = mdt_feature_atom_add(mlib, "Atom Z-coordinate", MOD_MDTC_NONE,
                               pos2, getbin, NULL, NULL);
  mdt_feature_add_needed_file(mlib, ifeat, MOD_MDTF_STRUCTURE);
  return ifeat;
}
