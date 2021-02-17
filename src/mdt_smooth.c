/** \file mdt_smooth.c     Functions to smooth MDTs.
 *
 *             Part of MDT, Copyright(c) 1989-2021 Andrej Sali
 */

#include <stdlib.h>
#include "modeller.h"
#include "mdt.h"
#include "util.h"

/** Smooth a histogram or the 2D plot with a uniform prior.
    Return TRUE on success. */
gboolean mdt_smooth(const struct mdt *mdtin, struct mdt *mdtout,
                    int dimensions, float weight, GError **err)
{
  static const float divisor = 1e-15;
  static const char *routine = "mdt_smooth";
  int nbins, nbinx, nbiny, *indf;

  if (!get_binx_biny(dimensions, &mdtin->base, routine, &nbinx, &nbiny, err)) {
    return FALSE;
  }
  nbins = nbinx * nbiny;

  mdt_copy(mdtin, mdtout, mdtin->base.bin_type);

  indf = mdt_start_indices(&mdtin->base);

  do {
    int i1, i2, i;
    float norm, w1, w2, wunifp;

    i1 = indmdt(indf, &mdtin->base);
    i2 = i1 + nbins;

    norm = 0.;
    for (i = i1; i < i2; i++) {
      norm += mod_mdt_bin_get(&mdtin->base, i);
    }

/*  The final distribution is: P = w1 * P(uniform) + w2 * P(data)
    where P_i(data) = F_i(data) / sum_i F_i(data) */

/*  weight is a constant which determines the number of points in
    F(data) at which the two input distributions have the same
    influence on the smoothed distributions. */
    weights(weight, nbins, norm, &w1, &w2);

/*  the weighted uniform p(i) */
    wunifp = w1 / (float)nbins;
/*  the weight for experimental frequencies */
    w2 = (norm > divisor ? w2 / norm : 0.);

    for (i = i1; i < i2; i++) {
      mod_mdt_bin_set(&mdtout->base, i,
                      wunifp + w2 * mod_mdt_bin_get(&mdtin->base, i));
    }

/* roll the indices of the "constant" features one forward: */
  } while (roll_ind_mdt(indf, &mdtin->base, mdtin->base.nfeat - dimensions));

  mdtout->pdf = TRUE;
  free(indf);
  return TRUE;
}
