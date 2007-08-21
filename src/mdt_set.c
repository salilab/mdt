/** \file mdt_set.c        Functions to set elements in MDTs.
 *
 *             Part of MDT, Copyright(c) 1989-2007 Andrej Sali
 */

#include <stdlib.h>
#include "modeller.h"
#include "mdt.h"
#include "util.h"

/** Set an element in an MDT. Return TRUE on success. */
gboolean mdt_set(struct mod_mdt *mdt, const int indices[], int n_indices,
                 double val, GError **err)
{
  const static char *routine = "mdt_set";
  int i, *indf, indx;
  if (n_indices != mdt->nfeat) {
    g_set_error(err, MDT_ERROR, MDT_ERROR_VALUE,
                "%s: Number of indices (%d) must match dimension of MDT (%d)",
                routine, n_indices, mdt->nfeat);
    return FALSE;
  }
  indf = g_malloc(sizeof(int) * n_indices);
  for (i = 0; i < n_indices; i++) {
    /* count negative indices from the end of the feature, as in Python */
    if (indices[i] < 0) {
      indf[i] = mdt->features[i].iend + 1 + indices[i];
    } else {
      indf[i] = indices[i] + 1;
    }
  }
  indx = indmdt(indf, mdt);
  g_free(indf);
  if (indx < 0 || indx >= mdt->nelems) {
    g_set_error(err, MDT_ERROR, MDT_ERROR_INDEX,
                "%s: Index %d out of range %d to %d", routine, indx, 0,
                mdt->nelems);
    return FALSE;
  } else {
    mdt->bin[indx] = val;
    return TRUE;
  }
}
