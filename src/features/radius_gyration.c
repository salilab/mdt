/** \file radius_gyration.c     Protein radius of gyration feature.
 *
 *             Part of MDT, Copyright(c) 1989-2008 Andrej Sali
 */

#include "modeller.h"
#include "../mdt_index.h"
#include "../mdt_feature.h"
#include "../mdt_property.h"

static int getbin(const struct mod_alignment *aln, int protein,
                  struct mdt_properties *prop, void *data,
                  const struct mod_mdt_libfeature *feat,
                  const struct mdt_library *mlib,
                  const struct mod_libraries *libs, GError **err)
{
  float f = property_radius_gyration(aln, protein, prop);
  return iclsbin(f, feat);
}

int mdt_feature_radius_of_gyration(struct mdt_library *mlib, int protein,
                                   GError **err)
{
  int ifeat;
  ifeat = mdt_feature_protein_add(mlib, "Radius of gyration", MOD_MDTC_NONE,
                                  protein, getbin, NULL, err);
  mdt_feature_add_needed_file(mlib, ifeat, MOD_MDTF_STRUCTURE);
  return ifeat;
}
