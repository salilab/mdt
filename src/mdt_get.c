/** \file mdt_get.c        Functions to get elements from MDTs.
 *
 *             Part of MDT, Copyright(c) 1989-2007 Andrej Sali
 */

#include <stdlib.h>
#include "modeller.h"
#include "mdt.h"
#include "util.h"

/** Get an element from an MDT. */
double mdt_get(const struct mdt_type *mdt, const int indices[], int n_indices,
               GError **err)
{
  const static char *routine = "mdt_get";
  int i, *indf, indx;
  if (n_indices != mdt->nfeat) {
    g_set_error(err, MDT_ERROR, MDT_ERROR_VALUE,
                "%s: Number of indices (%d) must match dimension of MDT (%d)",
                routine, n_indices, mdt->nfeat);
    return 0.0;
  }
  indf = g_malloc(sizeof(int) * n_indices);
  for (i = 0; i < n_indices; i++) {
    indf[i] = indices[i] + 1;
  }
  indx = indmdt(indf, mdt);
  g_free(indf);
  if (indx < 0 || indx >= mdt->nelems) {
    g_set_error(err, MDT_ERROR, MDT_ERROR_INDEX,
                "%s: Index %d out of range %d to %d", routine, indx, 0,
                mdt->nelems);
    return 0.0;
  } else {
    return mdt->bin[indx];
  }
}
