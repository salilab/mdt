/** \file mdt_linear_transform.c  Functions to transform MDTs.
 *
 *             Part of MDT, Copyright(c) 1989-2016 Andrej Sali
 */

#include "modeller.h"
#include "mdt.h"

/** Transform an MDT with a linear function */
void mdt_linear_transform(struct mod_mdt *mdt, float offset, float multiplier)
{
  int i;

  mod_lognote("transform_mdt_> parameters: %10.5f %10.5f\n"
              "                y = a + b*y", offset, multiplier);

  for (i = 0; i < mdt->nelems; i++) {
    mod_mdt_bin_set(mdt, i, offset + multiplier * mod_mdt_bin_get(mdt, i));
  }
}
