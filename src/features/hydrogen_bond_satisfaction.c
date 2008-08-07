/** \file hydrogen_bond_satisfaction.c
 *  \brief Protein hydrogen bond satisfaction feature.
 *
 *             Part of MDT, Copyright(c) 1989-2008 Andrej Sali
 */

#include "modeller.h"
#include "../mdt_index.h"
#include "../mdt_feature.h"
#include "../mdt_property.h"
#include "../mdt_hydrogen_bonds.h"

static int getbin(const struct mod_alignment *aln, int protein,
                  struct mdt_properties *prop, void *data,
                  const struct mod_mdt_libfeature *feat,
                  const struct mdt_library *mlib,
                  const struct mod_libraries *libs, GError **err)
{
  float fprop;
  if (property_hbpot(aln, protein, prop, mlib, libs, &fprop, err)) {
    return iclsbin(fprop, feat);
  } else {
    return -1;
  }
}

int mdt_feature_hydrogen_bond_satisfaction(struct mdt_library *mlib,
                                           int protein, GError **err)
{
  return mdt_feature_protein_add(mlib, "Protein hydrogen bond satisfaction",
                                 MOD_MDTC_NONE, protein, getbin, NULL, err);
}
