/** \file mdt_make.c       Functions to make new MDTs.
 *
 *             Part of MDT, Copyright(c) 1989-2011 Andrej Sali
 */

#include <glib.h>
#include <string.h>
#include "modeller.h"
#include "mdt.h"
#include "util.h"

/** Clear the MDT array, and set feature types. Return TRUE on success. */
gboolean mdt_make(struct mdt *mdt, const struct mdt_library *mlib,
                  const int features[], int n_features, GError **err)
{
  int i, nelems;

  nelems = 1;
  for (i = 0; i < n_features; i++) {
    int ifeat = features[i];
    if (!check_feature_type(ifeat, mlib, err)) {
      return FALSE;
    }
    nelems *= mlib->base.features[ifeat - 1].nbins;
  }
  mod_mdt_nelems_set(&mdt->base, nelems);
  mod_mdt_nfeat_set(&mdt->base, n_features);
  memset(mdt->base.bindata, 0, mod_mdt_bin_get_size(&mdt->base) * nelems);
  mdt->nalns = mdt->n_protein_pairs = mdt->n_proteins = 0;
  mdt->sample_size = 0.;

  for (i = 0; i < n_features; i++) {
    int ifeat = features[i];
    struct mod_mdt_feature *feat = &mdt->base.features[i];
    feat->ifeat = ifeat;
    feat->istart = 1;
    feat->iend = mlib->base.features[ifeat - 1].nbins;
  }

  return mdt_setup(mdt, mlib, err);
}
