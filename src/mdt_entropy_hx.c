/** \file mdt_entropy_hx.c  Functions to calculate entropy of dependent variable
 *
 *             Part of MDT, Copyright(c) 1989-2007 Andrej Sali
 */

#include <stdlib.h>
#include "modeller.h"
#include "mdt.h"
#include "util.h"

/** Get the entropy of the dependent variable. */
float mdt_entropy_hx(const struct mdt_type *mdt, GError **err)
{
  int i, nbinx;
  double summdt, sumfrq, *frq, hx;
  static const double small = 1.e-8;
  const char *routine = "mdt_entropy_hx";

  nbinx = f_int1_get(&mdt->nbins, mdt->nfeat - 1);

  /* the number of points in mdt */
  summdt = get_sum(f_double1_pt(&mdt->bin), mdt->nelems);

  if (summdt < small) {
    g_set_error(err, MDT_ERROR, MDT_ERROR_FAILED,
                "%s: MDT is empty; sum over all elements = %.4g", routine,
                summdt);
    return 0;
  }

  /* get pdf p(x) irrespective of the values of the independent variables */
  frq = g_malloc(sizeof(double) * nbinx);
  getfrq(mdt, NULL, 0, NULL, nbinx, frq);

  /* get its entropy: */
  sumfrq = get_sum(frq, nbinx);
  if (sumfrq < small) {
    sumfrq = small;
  }
  for (i = 0; i < nbinx; i++) {
    frq[i] = frq[i] / sumfrq;
  }

  hx = entrp1(frq, nbinx);
  if (hx < small) {
    modlogwarning(routine, "Entropy too small for division; changed to %.4g",
                  small);
    hx = small;
  }

  modlognote("entropy_hx_mdt_> h(x) = %.4g", hx);
  g_free(frq);
  return hx;
}
