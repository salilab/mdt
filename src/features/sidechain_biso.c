/** \file sidechain_biso.c Average sidechain Biso feature.
 *
 *             Part of MDT, Copyright(c) 1989-2011 Andrej Sali
 */

#include "modeller.h"
#include "../mdt_index.h"
#include "../mdt_feature.h"
#include "../mdt_all_features.h"
#include "../mdt_property.h"

static int getbin(const struct mod_alignment *aln, int protein, int residue,
                  struct mdt_properties *prop, void *data,
                  const struct mod_mdt_libfeature *feat,
                  const struct mdt_library *mlib,
                  const struct mod_libraries *libs, GError **err)
{
  const float *table;
  table = property_sidechain_biso(aln, protein, prop);
  if (table[residue] == 0.) {
    /* Biso of zero counts as undefined */
    return feat->nbins;
  } else {
    return feat_to_bin(table[residue], feat);
  }
}

int mdt_feature_sidechain_biso(struct mdt_library *mlib, int protein,
                               int delta, int align_delta, gboolean pos2,
                               GError **err)
{
  int ifeat;
  ifeat = mdt_feature_residue_add(mlib, "Average sidechain Biso", MOD_MDTC_NONE,
                                  protein, delta, align_delta, pos2, -1, getbin,
                                  NULL, NULL, err);
  if (ifeat < 0) {
    return ifeat;
  }
  mdt_feature_add_needed_file(mlib, ifeat, MOD_MDTF_STRUCTURE);
  return ifeat;
}
