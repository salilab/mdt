/** \file radius_gyration.c     Protein radius of gyration feature.
 *
 *             Part of MDT, Copyright(c) 1989-2015 Andrej Sali
 */

#include "modeller.h"
#include "../mdt_index.h"
#include "../mdt_feature.h"
#include "../mdt_all_features.h"
#include "../mdt_property.h"

static int getbin(const struct mod_alignment *aln, int protein,
                  struct mdt_properties *prop,
                  const struct mdt_feature *feat,
                  const struct mdt_library *mlib,
                  const struct mod_libraries *libs, GError **err)
{
  float f = property_radius_gyration(aln, protein, prop);
  return feat_to_bin(f, feat);
}

int mdt_feature_radius_of_gyration(struct mdt_library *mlib, int protein,
                                   GError **err)
{
  int ifeat;
  ifeat = mdt_feature_protein_add(mlib, "Radius of gyration", MOD_MDTC_NONE,
                                  protein, getbin, NULL, NULL, err);
  if (ifeat >= 0) {
    mdt_feature_add_needed_file(mlib, ifeat, MOD_MDTF_STRUCTURE);
  }
  return ifeat;
}
