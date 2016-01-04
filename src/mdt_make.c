/** \file mdt_make.c       Functions to make new MDTs.
 *
 *             Part of MDT, Copyright(c) 1989-2016 Andrej Sali
 */

#include <glib.h>
#include <string.h>
#include "modeller.h"
#include "mdt.h"
#include "util.h"

/** Clear the MDT array, and set feature types. Return TRUE on success. */
gboolean mdt_make(struct mdt *mdt, const struct mdt_library *mlib,
                  const int features[], int n_features,
                  const int shape[], int n_shape, GError **err)
{
  const char *routine = "mdt_reshape";
  int i, nelems;

  if (n_shape != n_features && n_shape != 0) {
    g_set_error(err, MDT_ERROR, MDT_ERROR_VALUE,
                "%s: shape dimension (%d) must be zero or match the number of "
                "features (%d)", routine, n_shape, n_features);
    return FALSE;
  }

  nelems = 1;
  for (i = 0; i < n_features; i++) {
    int ifeat = features[i];
    int nbins;

    if (!check_feature_type(ifeat, mlib, err)) {
      return FALSE;
    }

    nbins = mlib->base.features[ifeat - 1].nbins;
    if (n_shape > 0) {
      if (shape[i] <= -nbins || shape[i] > nbins) {
        g_set_error(err, MDT_ERROR, MDT_ERROR_VALUE,
                    "%s: shape[%d] (%d) is out of range - it must be between "
                    "%d and %d", routine, i, shape[i], 1 - nbins, nbins);
        return FALSE;
      }
      if (shape[i] <= 0) {
        nbins += shape[i];
      } else {
        nbins = shape[i];
      }
    }
    nelems *= nbins;
  }
  mod_mdt_nelems_set(&mdt->base, nelems);
  mod_mdt_nfeat_set(&mdt->base, n_features);
  memset(mdt->base.bindata, 0, mod_mdt_bin_get_size(&mdt->base) * nelems);
  mdt->nalns = mdt->n_protein_pairs = mdt->n_proteins = 0;
  mdt->sample_size = 0.;

  for (i = 0; i < n_features; i++) {
    int ifeat = features[i];
    int nbins = mlib->base.features[ifeat - 1].nbins;
    struct mod_mdt_feature *feat = &mdt->base.features[i];
    if (n_shape > 0) {
      if (shape[i] <= 0) {
        nbins += shape[i];
      } else {
        nbins = shape[i];
      }
    }
    feat->ifeat = ifeat;
    feat->istart = 1;
    feat->iend = nbins;
  }

  return mdt_setup(mdt, mlib, err);
}
