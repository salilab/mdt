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
  double *bin;

  modlognote("transform_mdt_> parameters: %10.5f %10.5f\n"
             "                y = Ln[a + b*y]", offset, multiplier);

  bin = f_double1_pt(&mdt->bin);

  for (i = 0; i < mdt->nelems; i++) {
    double arg = offset + multiplier * bin[i];
    if (arg < 1e-10) {
      bin[i] = undefined;
    } else {
      bin[i] = log(arg);
    }
  }
}
