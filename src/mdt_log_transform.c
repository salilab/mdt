/** \file mdt_log_transform.c  Functions to transform MDTs.
 *
 *             Part of MDT, Copyright(c) 1989-2011 Andrej Sali
 */

#include <stdlib.h>
#include <math.h>
#include "modeller.h"
#include "mdt.h"

/** Transform an MDT with a log function */
void mdt_log_transform(struct mod_mdt *mdt, float offset, float multiplier,
                       float undefined)
{
  int i;

  mod_lognote("transform_mdt_> parameters: %10.5f %10.5f\n"
              "                y = Ln[a + b*y]", offset, multiplier);

  for (i = 0; i < mdt->nelems; i++) {
    double arg = offset + multiplier * mod_mdt_bin_get(mdt, i);
    if (arg < 1e-10) {
      mod_mdt_bin_set(mdt, i, undefined);
    } else {
      mod_mdt_bin_set(mdt, i, log(arg));
    }
  }
}
