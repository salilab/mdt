/** \file residue_type.c  Residue type feature.
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
  struct mod_sequence *seq = mod_alignment_sequence_get(aln, protein);
  int restyp = mod_int1_get(&seq->irestyp, residue);
  if (restyp >= 21 || restyp <= 0) {
    mod_logwarning("residue_type", "Non-standard residue type (%d) at "
                   "position %d in protein %d - binning as undefined",
                   restyp, residue, protein);
    return feat->nbins;
  } else {
    return restyp;
  }
}

int mdt_feature_residue_type(struct mdt_library *mlib, int protein,
                             int delta, int align_delta, gboolean pos2,
                             const struct mod_libraries *libs, GError **err)
{
  int ifeat, i;
  struct mod_mdt_libfeature *feat;
  const static int nrestyp = 21;
  ifeat = mdt_feature_residue_add(mlib, "Residue type", MOD_MDTC_NONE,
                                  protein, delta, align_delta, pos2,
                                  libs->igaptyp, getbin, NULL, NULL, err);
  if (ifeat < 0) {
    return ifeat;
  }

  /* Set bins */
  feat = &mlib->base.features[ifeat - 1];
  mod_mdt_libfeature_nbins_set(feat, nrestyp + 1);
  for (i = 0; i < nrestyp; ++i) {
    g_free(feat->bins[i].symbol);
    feat->bins[i].symbol = mod_residue_name_from_type(i + 1, libs);
    feat->bins[i].rang1 = i;
    feat->bins[i].rang2 = i + 1;
  }
  g_free(feat->bins[nrestyp].symbol);
  feat->bins[nrestyp].symbol = g_strdup("u");
  feat->bins[nrestyp].rang1 = nrestyp;
  feat->bins[nrestyp].rang2 = nrestyp + 1;
  return ifeat;
}
