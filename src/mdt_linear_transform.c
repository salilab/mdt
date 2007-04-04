/** \file mdt_linear_transform.c  Functions to transform MDTs.
 *
 *             Part of MDT, Copyright(c) 1989-2007 Andrej Sali
 */

#include "modeller.h"
#include "mdt.h"

/** Transform an MDT with a linear function */
void mdt_linear_transform(struct mdt_type *mdt, float offset, float multiplier)
{
  int i;
  double *bin;

  modlognote("transform_mdt_> parameters: %10.5f %10.5f\n"
             "                y = a + b*y", offset, multiplier);

  bin = f_double1_pt(&mdt->bin);

  for (i = 0; i < mdt->nelems; i++) {
    bin[i] = offset + multiplier * bin[i];
  }
}
