/** \file mdt_exp_transform.c  Functions to transform MDTs.
 *
 *             Part of MDT, Copyright(c) 1989-2007 Andrej Sali
 */

#include <math.h>
#include "modeller.h"
#include "mdt.h"

/** Transform an MDT with an exp function */
void mdt_exp_transform(struct mdt_type *mdt, float offset, float expoffset,
                       float multiplier, float power)
{
  int i;
  double *bin;

  modlognote("transform_mdt_> parameters: %10.5f %10.5f %10.5f %10.5f\n"
             "                y = a + exp[b + c*y^d]", offset, expoffset,
             multiplier, power);

  bin = f_double1_pt(&mdt->bin);

  for (i = 0; i < mdt->nelems; i++) {
    bin[i] = offset + exp(expoffset + multiplier * pow(bin[i], power));
  }
}
