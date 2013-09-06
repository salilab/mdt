/** \file residue_group.c  Residue group feature.
 *
 *             Part of MDT, Copyright(c) 1989-2013 Andrej Sali
 */

#include "modeller.h"
#include "../mdt_error.h"
#include "../mdt_index.h"
#include "../mdt_feature.h"
#include "../mdt_all_features.h"

static int getbin(const struct mod_alignment *aln, int protein, int residue,
                  struct mdt_properties *prop,
                  const struct mdt_feature *feat,
                  const struct mdt_library *mlib,
                  const struct mod_libraries *libs, GError **err)
{
  int residue_grouping = GPOINTER_TO_INT(feat->data);
  struct mod_sequence *s = mod_alignment_sequence_get(aln, protein);
  int restyp = mod_int1_get(&s->irestyp, residue);
  return mod_residue_group_from_type(restyp, residue_grouping, libs);
}

int mdt_feature_residue_group(struct mdt_library *mlib, int protein, int delta,
                              int align_delta, gboolean pos2,
                              int residue_grouping,
                              const struct mod_libraries *libs, GError **err)
{
  int ifeat, i, nrescls;
  struct mod_mdt_libfeature *feat;
  char *name;

  if (residue_grouping < 0 || residue_grouping >= libs->nresgrp) {
    g_set_error(err, MDT_ERROR, MDT_ERROR_VALUE,
                "residue_grouping (%d) out of range; must be between 0 and %d",
                residue_grouping, libs->nresgrp - 1);
    return -1;
  }

  name = g_strdup_printf("Residue group type (residue_grouping=%d)",
                         residue_grouping);
  ifeat = mdt_feature_residue_add(mlib, name, MOD_MDTC_NONE, protein, delta,
                                  align_delta, pos2, -1, getbin,
                                  GINT_TO_POINTER(residue_grouping), NULL, err);
  g_free(name);
  if (ifeat < 0) {
    return ifeat;
  }

  /* Set bins */
  nrescls = libs->nrescls[residue_grouping];
  feat = &mlib->base.features[ifeat - 1];
  mod_mdt_libfeature_nbins_set(feat, nrescls + 1);
  for (i = 0; i < nrescls; ++i) {
    g_free(feat->bins[i].symbol);
    feat->bins[i].symbol = g_strdup_printf("%d", i);
    feat->bins[i].rang1 = i;
    feat->bins[i].rang2 = i + 1;
  }
  g_free(feat->bins[nrescls].symbol);
  feat->bins[nrescls].symbol = g_strdup("U");
  feat->bins[nrescls].rang1 = nrescls;
  feat->bins[nrescls].rang2 = nrescls + 1;
  return ifeat;
}
