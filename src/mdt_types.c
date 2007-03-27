/** \file mdt_types.c  Functions to handle MDT types.
 *
 *             Part of MDT, Copyright(c) 1989-2007 Andrej Sali
 */

#include <glib.h>
#include "modeller.h"
#include "mdt_types.h"

/** Make a new mdt_library structure */
struct mdt_library *mdt_library_new(void)
{
  struct mdt_library *mlib;
  mlib = g_malloc(sizeof(struct mdt_library));
  mod_mdt_library_init(&mlib->base);
  mlib->deltai = mlib->deltaj = 1;
  mlib->deltai_ali = mlib->deltaj_ali = FALSE;
  return mlib;
}

/** Free an mdt_library structure */
void mdt_library_free(struct mdt_library *mlib)
{
  mod_mdt_library_dealloc(&mlib->base);
  g_free(mlib);
}
