/** \file mdt_close.c   Functions to close MDTs, e.g. for splines.
 *
 *             Part of MDT, Copyright(c) 1989-2007 Andrej Sali
 */

#include <stdlib.h>
#include "modeller.h"
#include "mdt.h"
#include "util.h"

/** Close a 1D array */
static void close_1d(double *bin, int nbins)
{
  bin[0] = 0.5 * (bin[0] + bin[nbins - 1]);
  bin[nbins - 1] = bin[0];
}


/** Close a 2D array. Elements are addressed as bin[y*nbinx + x] */
static void close_2d(double *bin, int nbinx, int nbiny)
{
  int x, y;
  double mean;

  /* The vertical edges must be the same: */
  for (x = 0; x < nbinx; x++) {
    mean = 0.5 * (bin[0 * nbinx + x] + bin[(nbiny - 1) * nbinx + x]);
    bin[0 * nbinx + x] = bin[(nbiny - 1) * nbinx + x] = mean;
  }

  /* The horizontal edges must be the same: */
  for (y = 0; y < nbiny; y++) {
    mean = 0.5 * (bin[y * nbinx + 0] + bin[y * nbinx + nbinx - 1]);
    bin[y * nbinx + 0] = bin[y * nbinx + nbinx - 1] = mean;
  }

  /* The four corners must be the same: */
  mean = 0.25 * (bin[0 * nbinx + 0]
                 + bin[(nbiny - 1) * nbinx + 0]
                 + bin[0 * nbinx + nbinx - 1]
                 + bin[(nbiny - 1) * nbinx + nbinx - 1]);
  bin[0 * nbinx + 0] = bin[(nbiny - 1) * nbinx + 0]
      = bin[0 * nbinx + nbinx - 1]
      = bin[(nbiny - 1) * nbinx + nbinx - 1] = mean;
}


/** Close an MDT so that it is useful for creating periodic splines.
    Return TRUE on success. */
gboolean mdt_close(struct mdt_type *mdt, int dimensions, GError **err)
{
  static const char *routine = "mdt_close";
  int nbins, nbinx, nbiny, *indf;
  double *bin;

  if (!get_binx_biny(dimensions, mdt, routine, &nbinx, &nbiny, err)) {
    return FALSE;
  }
  nbins = nbinx * nbiny;

  modlognote("transform_mdt_> close the ends of a spline");

  indf = mdt_start_indices(mdt);
  bin = f_double1_pt(&mdt->bin);

  do {
    int i1 = indmdt(indf, mdt);
    if (dimensions == 1) {
      close_1d(&bin[i1], nbins);
    } else {
      close_2d(&bin[i1], nbinx, nbiny);
    }
/* roll the indices of the "constant" features one forward: */
  } while (roll_ind(indf, f_int1_pt(&mdt->istart), f_int1_pt(&mdt->iend),
                    mdt->nfeat - dimensions));

  free(indf);
  return TRUE;
}
