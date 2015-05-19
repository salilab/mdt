/** \file atom_accessibility.c  Atom accessibility feature.
 *
 *             Part of MDT, Copyright(c) 1989-2015 Andrej Sali
 */

#include "modeller.h"
#include "../mdt_index.h"
#include "../mdt_feature.h"
#include "../mdt_all_features.h"
#include "../geometry.h"

static int getbin(const struct mod_alignment *aln, int protein, int atom,
                  struct mdt_properties *prop,
                  const struct mdt_feature *feat,
                  const struct mdt_library *mlib,
                  const struct mod_libraries *libs, GError **err)
{
  struct mod_structure *s = mod_alignment_structure_get(aln, protein);
  float f = mod_float1_get(&s->cd.x, atom);
  if (coordinate_undefined(f)) {
    return mdt_feature_undefined_bin_get(feat);
  } else {
    f = mod_float1_get(&s->cd.atmacc, atom);
    return feat_to_bin(f, feat);
  }
}

int mdt_feature_atom_accessibility(struct mdt_library *mlib, gboolean pos2)
{
  int ifeat;
  ifeat = mdt_feature_atom_add(mlib, "Atom accessibility", MOD_MDTC_NONE,
                               pos2, getbin, NULL, NULL);
  mdt_feature_add_needed_file(mlib, ifeat, MOD_MDTF_PSA);
  return ifeat;
}
