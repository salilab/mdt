/** \file mdt_inverse_transform.c  Functions to transform MDTs.
 *
 *             Part of MDT, Copyright(c) 1989-2007 Andrej Sali
 */

#include <stdlib.h>
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
    if (abs(mdt->bin[i]) < divisor) {
      mdt->bin[i] = undefined;
    } else {
      mdt->bin[i] = offset + multiplier / mdt->bin[i];
    }
  }
}
