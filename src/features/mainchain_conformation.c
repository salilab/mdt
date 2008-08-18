/** \file mainchain_conformation.c  Ramachandran mainchain conformation feature.
 *
 *             Part of MDT, Copyright(c) 1989-2008 Andrej Sali
 */

#include "modeller.h"
#include "../mdt_index.h"
#include "../mdt_feature.h"
#include "../mdt_all_features.h"

static int getbin(const struct mod_alignment *aln, int protein, int residue,
                  struct mdt_properties *prop, void *data,
                  const struct mod_mdt_libfeature *feat,
                  const struct mdt_library *mlib,
                  const struct mod_libraries *libs, GError **err)
{
  struct mod_structure *s = mod_alignment_structure_get(aln, protein);
  return mod_int1_get(&s->imnchw, residue);
}

int mdt_feature_mainchain_conformation(struct mdt_library *mlib, int protein,
                                       int delta, int align_delta,
                                       gboolean pos2,
                                       const struct mod_libraries *libs,
                                       GError **err)
{
  int ifeat, i;
  struct mod_mdt_libfeature *feat;
  ifeat = mdt_feature_residue_add(mlib, "Mainchain conformation (Ramachandran)",
                                  MOD_MDTC_MNRAMA, protein, delta, align_delta,
                                  pos2, -1, getbin, NULL, NULL, err);
  if (ifeat < 0) {
    return ifeat;
  }
  mdt_feature_add_needed_file(mlib, ifeat, MOD_MDTF_STRUCTURE);
  mdt_feature_add_needed_file(mlib, ifeat, MOD_MDTF_DIHEDRALS);

  /* Set bins */
  feat = &mlib->base.features[ifeat - 1];
  mod_mdt_libfeature_nbins_set(feat, libs->nmnch + 1);
  for (i = 0; i < libs->nmnch; ++i) {
    g_free(feat->bins[i].symbol);
    feat->bins[i].symbol = g_strdup_printf("%c", libs->wilmot.mnchcls[i]);
    feat->bins[i].rang1 = i;
    feat->bins[i].rang2 = i + 1;
  }
  g_free(feat->bins[libs->nmnch].symbol);
  feat->bins[libs->nmnch].symbol = g_strdup("U");
  feat->bins[libs->nmnch].rang1 = libs->nmnch;
  feat->bins[libs->nmnch].rang2 = libs->nmnch + 1;
  return ifeat;
}
