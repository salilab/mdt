/** \file mdt_log_transform.c  Functions to transform MDTs.
 *
 *             Part of MDT, Copyright(c) 1989-2007 Andrej Sali
 */

#include <stdlib.h>
#include <math.h>
#include "modeller.h"
#include "mdt.h"

/** Transform an MDT with a log function */
void mdt_log_transform(struct mdt_type *mdt, float offset, float multiplier,
                       float undefined)
{
  int i;

  modlognote("transform_mdt_> parameters: %10.5f %10.5f\n"
             "                y = Ln[a + b*y]", offset, multiplier);

  for (i = 0; i < mdt->nelems; i++) {
    double arg = offset + multiplier * mdt->bin[i];
    if (arg < 1e-10) {
      mdt->bin[i] = undefined;
    } else {
      mdt->bin[i] = log(arg);
    }
  }
}
