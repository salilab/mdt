/** \file fractional_atom_accessibility.c Fractional atom accessibility feature.
 *
 *             Part of MDT, Copyright(c) 1989-2008 Andrej Sali
 */

#include "modeller.h"
#include "../mdt_index.h"
#include "../mdt_feature.h"
#include "../mdt_property.h"

static int getbin(const struct mod_alignment *aln, int protein, int atom,
                  struct mdt_properties *prop, void *data,
                  const struct mod_mdt_libfeature *feat,
                  const struct mdt_library *mlib,
                  const struct mod_libraries *libs, GError **err)
{
  float *table;
  struct mod_structure *s = mod_alignment_structure_get(aln, protein);
  if (property_fatmacc(aln, protein, prop, feat, libs, &table, err)) {
    return ftable(table, s->cd.natm, atom, feat);
  } else {
    return -1;
  }
}

int mdt_feature_fractional_atom_accessibility(struct mdt_library *mlib,
                                              gboolean pos2)
{
  int ifeat;
  ifeat = mdt_feature_atom_add(mlib, "Fractional atom accessibility",
                               MOD_MDTC_NONE, pos2, getbin, NULL);
  mdt_feature_add_needed_file(mlib, ifeat, MOD_MDTF_STRUCTURE);
  mdt_feature_add_needed_file(mlib, ifeat, MOD_MDTF_PSA);
  return ifeat;
}
