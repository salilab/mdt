/** \file residue_accessibility.c  Residue accessibility feature.
 *
 *             Part of MDT, Copyright(c) 1989-2008 Andrej Sali
 */

#include "modeller.h"
#include "../mdt_index.h"
#include "../mdt_feature.h"

static int getbin(const struct mod_alignment *aln, int protein, int residue,
                  int alnpos, struct mdt_properties *prop, void *data,
                  const struct mod_mdt_libfeature *feat,
                  const struct mdt_library *mlib,
                  const struct mod_libraries *libs, GError **err)
{
  struct mod_structure *s = mod_alignment_structure_get(aln, protein);
  float f = mod_float1_get(&s->acc, residue);
  return iclsbin(f, feat);
}

int mdt_feature_residue_accessibility(struct mdt_library *mlib, int protein,
                                      int delta, gboolean pos2, GError **err)
{
  int ifeat;
  ifeat = mdt_feature_residue_add(mlib, "Residue accessibility", MOD_MDTC_NONE,
                                  protein, delta, pos2, -1, getbin, NULL, err);
  mdt_feature_add_needed_file(mlib, ifeat, MOD_MDTF_STRUCTURE);
  mdt_feature_add_needed_file(mlib, ifeat, MOD_MDTF_PSA);
  return ifeat;
}
