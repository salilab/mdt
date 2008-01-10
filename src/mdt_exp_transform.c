/** \file mdt_exp_transform.c  Functions to transform MDTs.
 *
 *             Part of MDT, Copyright(c) 1989-2008 Andrej Sali
 */

#include <math.h>
#include "modeller.h"
#include "mdt.h"

/** Transform an MDT with an exp function */
void mdt_exp_transform(struct mod_mdt *mdt, float offset, float expoffset,
                       float multiplier, float power)
{
  int i;

  mod_lognote("transform_mdt_> parameters: %10.5f %10.5f %10.5f %10.5f\n"
              "                y = a + exp[b + c*y^d]", offset, expoffset,
              multiplier, power);

  for (i = 0; i < mdt->nelems; i++) {
    mdt->bin[i] = offset + exp(expoffset
                               + multiplier * pow(mdt->bin[i], power));
  }
}
