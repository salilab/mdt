/** \file mdt_clear.c    Functions to set an MDT to zero.
 *
 *             Part of MDT, Copyright(c) 1989-2016 Andrej Sali
 */

#include <string.h>
#include "mdt.h"

/** Clear an MDT (set all elements to zero). */
void mdt_clear(struct mdt *mdt)
{
  memset(mdt->base.bindata, 0,
         mod_mdt_bin_get_size(&mdt->base) * mdt->base.nelems);
}
