/** \file mdt_offset_min.c   Functions to offset MDTs by the minimum value.
 *
 *             Part of MDT, Copyright(c) 1989-2013 Andrej Sali
 */

#include <stdlib.h>
#include "modeller.h"
#include "mdt.h"
#include "util.h"

/** Get the minimum value in the bin array */
static double get_bin_minval(const struct mod_mdt *mdt, int offset, int num)
{
  int i;
  double minval;
  if (num == 0) {
    return 0;
  }
  minval = mod_mdt_bin_get(mdt, offset);
  for (i = 1; i < num; i++) {
    double binval = mod_mdt_bin_get(mdt, offset + i);
    if (binval < minval) {
      minval = binval;
    }
  }
  return minval;
}

/** Offset an MDT by the minimum value. Return TRUE on success. */
gboolean mdt_offset_min(struct mod_mdt *mdt, int dimensions, GError **err)
{
  static const char *routine = "mdt_offset_min";
  int nbins, nbinx, nbiny, *indf;

  if (!get_binx_biny(dimensions, mdt, routine, &nbinx, &nbiny, err)) {
    return FALSE;
  }
  nbins = nbinx * nbiny;

  mod_lognote("transform_mdt_> parameters:\n"
              "                y = y - min(y)");

  indf = mdt_start_indices(mdt);

  do {
    int i1, i2, i;
    double minval;

    i1 = indmdt(indf, mdt);
    i2 = i1 + nbins;
    minval = get_bin_minval(mdt, i1, nbins);

    for (i = i1; i < i2; i++) {
      mod_mdt_bin_set(mdt, i, mod_mdt_bin_get(mdt, i) - minval);
    }

/* roll the indices of the "constant" features one forward: */
  } while (roll_ind_mdt(indf, mdt, mdt->nfeat - dimensions));

  free(indf);
  return TRUE;
}
