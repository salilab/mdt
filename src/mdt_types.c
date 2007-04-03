/** \file mdt_types.c  Functions to handle MDT types.
 *
 *             Part of MDT, Copyright(c) 1989-2007 Andrej Sali
 */

#include <glib.h>
#include "modeller.h"
#include "mdt_types.h"
#include "mdt_atom_classes.h"

/** Make a new mdt_library structure */
struct mdt_library *mdt_library_new(void)
{
  int i;
  struct mdt_library *mlib;
  mlib = g_malloc(sizeof(struct mdt_library));
  mod_mdt_library_init(&mlib->base);
  mlib->deltai = mlib->deltaj = 1;
  mlib->deltai_ali = mlib->deltaj_ali = FALSE;
  mlib->hbond_cutoff = 3.5;
  mlib->special_atoms = FALSE;
  for (i = 0; i < 4; i++) {
    mlib->atclass[i] = mdt_atom_class_list_new(i + 1);
  }
  mlib->hbond = mdt_atom_class_list_new(1);
  return mlib;
}

/** Free an mdt_library structure */
void mdt_library_free(struct mdt_library *mlib)
{
  int i;
  mod_mdt_library_dealloc(&mlib->base);
  for (i = 0; i < 4; i++) {
    mdt_atom_class_list_free(mlib->atclass[i]);
  }
  mdt_atom_class_list_free(mlib->hbond);
  g_free(mlib);
}
