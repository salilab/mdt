/** \file mdt_feature.c    Functions to act on MDT features.
 *
 *             Part of MDT, Copyright(c) 1989-2008 Andrej Sali
 */

#include <stdlib.h>
#include "modeller.h"
#include "mdt.h"
#include "mdt_feature.h"

/** Is the given feature type periodic? */
gboolean mdt_feature_is_periodic(int ifeat)
{
  switch (ifeat) {
  case 3:
  case 4:
  case 5:
  case 6:
  case 7:
  case 8:
  case 9:
  case 10:
  case 28:
  case 29:
  case 41:
  case 42:
  case 53:
  case 54:
  case 55:
  case 56:
  case 57:
  case 106:
  case 107:
  case 108:
  case 114:
    return TRUE;
  default:
    return FALSE;
  }
}

void mdt_feature_add_needed_file(struct mdt_library *mlib, int ifeat,
                                 mod_mdt_file filetype)
{
  struct mod_mdt_libfeature *feat = &mlib->base.features[ifeat - 1];
  feat->ndatfil++;
  feat->idatfil = g_realloc(feat->idatfil, sizeof(int) * feat->ndatfil);
  feat->idatfil[feat->ndatfil - 1] = filetype;
}

void mdt_feature_nbins_set(struct mdt_library *mlib, int ifeat,
                           int nbins)
{
  struct mod_mdt_libfeature *libfeat;
  libfeat = &mlib->base.features[ifeat - 1];

  /* Add one extra for the undefined bin */
  mod_mdt_libfeature_nbins_set(libfeat, nbins + 1);
  libfeat->bins[nbins].rang1 = 0.;
  libfeat->bins[nbins].rang2 = 1.;
  g_free(libfeat->bins[nbins].symbol);
  libfeat->bins[nbins].symbol = g_strdup("U");
}

void mdt_feature_bin_set(struct mdt_library *mlib, int ifeat, int bin,
                         float start, float end, const char *symbol)
{
  struct mod_mdt_libfeature *libfeat;
  libfeat = &mlib->base.features[ifeat - 1];

  libfeat->bins[bin].rang1 = start;
  libfeat->bins[bin].rang2 = end;
  g_free(libfeat->bins[bin].symbol);
  libfeat->bins[bin].symbol = g_strdup(symbol ? symbol : "");
}

int mdt_feature_protein_add(struct mdt_library *mlib, const char *name,
                            mod_mdt_calc precalc_type, int protein,
                            mdt_cb_feature_protein getbin, void *data,
                            GError **err)
{
  GString *fullname;
  struct mdt_feature *feat;
  int nfeat = mlib->base.nfeat + 1;

  if (protein < 0 || protein > 1) {
    g_set_error(err, MDT_ERROR, MDT_ERROR_VALUE,
                "Protein features can act only on protein 0 or protein 1; "
                "%d was given", protein);
    return -1;
  }

  mlib->features = g_array_set_size(mlib->features, nfeat);
  feat = &g_array_index(mlib->features, struct mdt_feature, nfeat - 1);
  feat->type = MDT_FEATURE_PROTEIN;
  feat->u.protein.protein = protein;
  feat->u.protein.getbin = getbin;
  feat->data = data;
  fullname = g_string_new(name);
  g_string_append_printf(fullname, " in protein %d", protein);
  mod_mdt_libfeature_register(&mlib->base, nfeat, fullname->str, precalc_type,
                              protein == 0 ? MOD_MDTP_A : MOD_MDTP_B,
                              MOD_MDTS_PROTEIN, FALSE, 0);
  g_string_free(fullname, TRUE);
  return nfeat;
}
