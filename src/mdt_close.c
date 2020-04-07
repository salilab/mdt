/** \file mdt_close.c   Functions to close MDTs, e.g. for splines.
 *
 *             Part of MDT, Copyright(c) 1989-2020 Andrej Sali
 */

#include <stdlib.h>
#include "modeller.h"
#include "mdt.h"
#include "util.h"

/** Close a 1D array */
static void close_1d(struct mod_mdt *mdt, int offset, int nbins)
{
  double binstart = mod_mdt_bin_get(mdt, offset);
  double binend = mod_mdt_bin_get(mdt, offset + nbins - 1);
  double binav = 0.5 * (binstart + binend);

  mod_mdt_bin_set(mdt, offset, binav);
  mod_mdt_bin_set(mdt, offset + nbins - 1, binav);
}


/** Close a 2D array. Elements are addressed as bin[offset + y*nbinx + x] */
static void close_2d(struct mod_mdt *mdt, int offset, int nbinx, int nbiny)
{
  int x, y;
  double mean;

  /* The vertical edges must be the same: */
  for (x = 0; x < nbinx; x++) {
    mean = 0.5 * (mod_mdt_bin_get(mdt, offset + 0 * nbinx + x)
                  + mod_mdt_bin_get(mdt, offset + (nbiny - 1) * nbinx + x));
    mod_mdt_bin_set(mdt, offset + 0 * nbinx + x, mean);
    mod_mdt_bin_set(mdt, offset + (nbiny - 1) * nbinx + x, mean);
  }

  /* The horizontal edges must be the same: */
  for (y = 0; y < nbiny; y++) {
    mean = 0.5 * (mod_mdt_bin_get(mdt, offset + y * nbinx + 0)
                  + mod_mdt_bin_get(mdt, offset + y * nbinx + nbinx - 1));
    mod_mdt_bin_set(mdt, offset + y * nbinx + 0, mean);
    mod_mdt_bin_set(mdt, offset + y * nbinx + nbinx - 1, mean);
  }

  /* The four corners must be the same: */
  mean = 0.25 * (mod_mdt_bin_get(mdt, offset + 0 * nbinx + 0)
                 + mod_mdt_bin_get(mdt, offset + (nbiny - 1) * nbinx + 0)
                 + mod_mdt_bin_get(mdt, offset + 0 * nbinx + nbinx - 1)
                 + mod_mdt_bin_get(mdt, offset + (nbiny - 1) * nbinx
                                        + nbinx - 1));
  mod_mdt_bin_set(mdt, offset + 0 * nbinx + 0, mean);
  mod_mdt_bin_set(mdt, offset + (nbiny - 1) * nbinx + 0, mean);
  mod_mdt_bin_set(mdt, offset + 0 * nbinx + nbinx - 1, mean);
  mod_mdt_bin_set(mdt, offset + (nbiny - 1) * nbinx + nbinx - 1, mean);
}


/** Close an MDT so that it is useful for creating periodic splines.
    Return TRUE on success. */
gboolean mdt_close(struct mod_mdt *mdt, int dimensions, GError **err)
{
  static const char *routine = "mdt_close";
  int nbins, nbinx, nbiny, *indf;

  if (!get_binx_biny(dimensions, mdt, routine, &nbinx, &nbiny, err)) {
    return FALSE;
  }
  nbins = nbinx * nbiny;

  mod_lognote("transform_mdt_> close the ends of a spline");

  indf = mdt_start_indices(mdt);

  do {
    int i1 = indmdt(indf, mdt);
    if (dimensions == 1) {
      close_1d(mdt, i1, nbins);
    } else {
      close_2d(mdt, i1, nbinx, nbiny);
    }
/* roll the indices of the "constant" features one forward: */
  } while (roll_ind_mdt(indf, mdt, mdt->nfeat - dimensions));

  g_free(indf);
  return TRUE;
}
