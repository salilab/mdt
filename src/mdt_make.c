/** \file mdt_make.c       Functions to make new MDTs.
 *
 *             Part of MDT, Copyright(c) 1989-2007 Andrej Sali
 */

#include <glib.h>
#include <string.h>
#include "modeller.h"
#include "mdt.h"
#include "util.h"

/** Clear the MDT array, and set feature types. Return TRUE on success. */
gboolean mdt_make(struct mdt_type *mdt, const struct mdt_library *mlib,
                  const int features[], int n_features, GError **err)
{
  const static char *routine = "mdt_make";
  int i, nelems, ierr;

  nelems = 1;
  for (i = 0; i < n_features; i++) {
    int ifeat = features[i];
    if (ifeat < 1 || ifeat > mlib->base.nfeat) {
      g_set_error(err, MDT_ERROR, MDT_ERROR_INDEX,
                  "%s: Feature type %d out of range 1-%d", routine, ifeat,
                  mlib->base.nfeat);
      return FALSE;
    }
    nelems *= mlib->base.features[ifeat - 1].nbins;
  }
  mdt_type_nelems_set(mdt, nelems);
  mdt_type_nfeat_set(mdt, n_features);
  memset(mdt->bin, 0, sizeof(double) * nelems);
  mdt->nalns = mdt->n_protein_pairs = mdt->n_proteins = 0;
  mdt->sample_size = 0.;

  for (i = 0; i < n_features; i++) {
    int ifeat = features[i];
    struct mdt_feature *feat = &mdt->features[i];
    feat->ifeat = ifeat;
    feat->istart = 1;
    feat->iend = mlib->base.features[ifeat - 1].nbins;
  }

  mdt_setup_check(mdt, &mlib->base, &ierr);
  if (ierr == 0) {
    return TRUE;
  } else {
    handle_modeller_error(err);
    return FALSE;
  }
}
