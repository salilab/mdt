/** \file radius_gyration.c     Protein radius of gyration feature.
 *
 *             Part of MDT, Copyright(c) 1989-2008 Andrej Sali
 */

#include "modeller.h"
#include "../mdt_feature.h"
#include "../mdt_property.h"

static float getfeat(const struct mod_alignment *aln, int protein,
                     struct mdt_properties *prop, void *data)
{
  return property_radius_gyration(aln, protein, prop);
}

int mdt_feature_radius_of_gyration(struct mdt_library *mlib, int protein,
                                   GError **err)
{
  return mdt_feature_protein_add(mlib, "Radius of gyration", MOD_MDTC_NONE,
                                 protein, getfeat, NULL, err);
}
