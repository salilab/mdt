/** \file mdt_smooth.c     Functions to smooth MDTs.
 *
 *             Part of MDT, Copyright(c) 1989-2006 Andrej Sali
 */

#include <stdlib.h>
#include "modeller.h"
#include "mdt.h"
#include "util.h"

/** Smooth a histogram or the 2D plot with a uniform prior */
void mdt_smooth(const struct mdt_type *mdtin, struct mdt_type *mdtout,
                int dimensions, float weight, int *ierr)
{
  static const float divisor = 1e-15;
  static const char *routine = "mdt_smooth";
  int nbins, nbinx, nbiny, *indf;
  double *in_bin, *out_bin;

  *ierr = 0;
  get_binx_biny(dimensions, mdtin, routine, &nbinx, &nbiny, ierr);
  if (*ierr != 0) {
    return;
  }
  nbins = nbinx * nbiny;

  copy_mdt(mdtin, mdtout);

  indf = mdt_start_indices(mdtin);
  in_bin = f_double1_pt(&mdtin->bin);
  out_bin = f_double1_pt(&mdtout->bin);

  do {
    int i1, i2, i;
    float norm, w1, w2, wunifp;

    i1 = indmdt(indf, mdtin);
    i2 = i1 + nbins;

    norm = 0.;
    for (i = i1; i < i2; i++) {
      norm += in_bin[i];
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
      out_bin[i] = wunifp + w2 * in_bin[i];
    }

/* roll the indices of the "constant" features one forward: */
  } while (roll_ind(indf, mdtin, mdtin->nfeat - dimensions));

  mdtout->pdf = 1;
  free(indf);
}
