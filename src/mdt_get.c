/** \file mdt_get.c        Functions to get elements from MDTs.
 *
 *             Part of MDT, Copyright(c) 1989-2006 Andrej Sali
 */

#include <stdlib.h>
#include "modeller.h"
#include "mdt.h"
#include "util.h"

/** Get an element from an MDT. */
double mdt_get(const struct mdt_type *mdt, const int indices[], int n_indices,
               int *ierr)
{
  const static char *routine = "mdt_get";
  int i, *indf, indx;
  *ierr = 0;
  if (n_indices != mdt->nfeat) {
    modlogerror(routine, ME_VALUE, "Number of indices (%d) must match "
                "dimension of MDT (%d)", n_indices, mdt->nfeat);
    *ierr = 1;
    return 0.0;
  }
  indf = malloc(sizeof(int) * n_indices);
  for (i = 0; i < n_indices; i++) {
    indf[i] = indices[i] + 1;
  }
  indx = indmdt(indf, mdt);
  free(indf);
  if (indx < 0 || indx > mdt->nelems) {
    modlogerror(routine, ME_INDEX, "Index %d out of range %d to %d", indx, 0,
                mdt->nelems);
    *ierr = 1;
    return 0.0;
  } else {
    double *bin = f_double1_pt(&mdt->bin);
    return bin[indx];
  }
}
