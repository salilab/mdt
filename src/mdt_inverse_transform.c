/** \file mdt_inverse_transform.c  Functions to transform MDTs.
 *
 *             Part of MDT, Copyright(c) 1989-2007 Andrej Sali
 */

#include <stdlib.h>
#include "modeller.h"
#include "mdt.h"

/** Transform an MDT with an inverse function */
void mdt_inverse_transform(struct mdt_type *mdt, float offset, float multiplier,
                           float undefined)
{
  static const float divisor = 1e-15;
  int i;
  double *bin;

  modlognote("transform_mdt_> parameters: %10.5f %10.5f\n"
             "                y = a + b/y", offset, multiplier);

  bin = f_double1_pt(&mdt->bin);

  for (i = 0; i < mdt->nelems; i++) {
    if (abs(bin[i]) < divisor) {
      bin[i] = undefined;
    } else {
      bin[i] = offset + multiplier / bin[i];
    }
  }
}
