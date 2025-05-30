/** \file atom_type.c  Modeller atom type feature.
 *
 *             Part of MDT, Copyright(c) 1989-2025 Andrej Sali
 */

#include "modeller.h"
#include "../mdt_index.h"
#include "../mdt_feature.h"
#include "../mdt_all_features.h"
#include "../mdt_property.h"
#include "../mdt_atom_classes.h"

static int getbin(const struct mod_alignment *aln, int protein, int atom,
                  struct mdt_properties *prop,
                  const struct mdt_feature *feat,
                  const struct mdt_library *mlib,
                  const struct mod_libraries *libs, GError **err)
{
  const int *binprop;
  binprop = property_iatta(aln, protein, prop, mlib, libs, err);
  if (binprop) {
    return binprop[atom];
  } else {
    return -1;
  }
}

int mdt_feature_atom_type(struct mdt_library *mlib, gboolean pos2)
{
  int ifeat;
  struct mod_mdt_libfeature *feat;
  ifeat = mdt_feature_atom_add(mlib, "Atom type", MOD_MDTC_NONE,
                               pos2, getbin, NULL, NULL);

  /* Set number of bins and their symbols */
  feat = &mlib->base.features[ifeat - 1];
  update_mdt_feat_atclass(feat, mlib->atclass[0]);
  return ifeat;
}
