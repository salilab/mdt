/** \file mdt_inverse_transform.c  Functions to transform MDTs.
 *
 *             Part of MDT, Copyright(c) 1989-2016 Andrej Sali
 */

#include <stdlib.h>
#include <math.h>
#include "modeller.h"
#include "mdt.h"

/** Transform an MDT with an inverse function */
void mdt_inverse_transform(struct mod_mdt *mdt, float offset,
                           float multiplier, float undefined)
{
  static const float divisor = 1e-15;
  int i;

  mod_lognote("transform_mdt_> parameters: %10.5f %10.5f\n"
              "                y = a + b/y", offset, multiplier);

  for (i = 0; i < mdt->nelems; i++) {
    double binval = mod_mdt_bin_get(mdt, i);
    if (fabs(binval) < divisor) {
      mod_mdt_bin_set(mdt, i, undefined);
    } else {
      mod_mdt_bin_set(mdt, i, offset + multiplier / binval);
    }
  }
}
