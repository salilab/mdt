/** \file fractional_atom_accessibility.c Fractional atom accessibility feature.
 *
 *             Part of MDT, Copyright(c) 1989-2025 Andrej Sali
 */

#include "modeller.h"
#include "../geometry.h"
#include "../mdt_index.h"
#include "../mdt_feature.h"
#include "../mdt_all_features.h"
#include "../mdt_property.h"

static int getbin(const struct mod_alignment *aln, int protein, int atom,
                  struct mdt_properties *prop,
                  const struct mdt_feature *feat,
                  const struct mdt_library *mlib,
                  const struct mod_libraries *libs, GError **err)
{
  const float *table;
  struct mod_structure *s = mod_alignment_structure_get(aln, protein);
  float f = mod_float1_get(&s->cd.x, atom);
  if (coordinate_undefined(f)) {
    return mdt_feature_undefined_bin_get(feat);
  } else if ((table = property_fatmacc(aln, protein, prop, libs, err))) {
    return feat_to_bin(table[atom], feat);
  } else {
    return -1;
  }
}

int mdt_feature_fractional_atom_accessibility(struct mdt_library *mlib,
                                              gboolean pos2)
{
  int ifeat;
  ifeat = mdt_feature_atom_add(mlib, "Fractional atom accessibility",
                               MOD_MDTC_NONE, pos2, getbin, NULL, NULL);
  mdt_feature_add_needed_file(mlib, ifeat, MOD_MDTF_PSA);
  return ifeat;
}
