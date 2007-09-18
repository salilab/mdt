/** \file mdt_read.c       Functions to read MDTs from text files.
 *
 *             Part of MDT, Copyright(c) 1989-2007 Andrej Sali
 */

#include <glib.h>
#include "modeller.h"
#include "mdt.h"
#include "util.h"

/** Read in an MDT in text format. Return TRUE on success. */
gboolean mdt_read(struct mdt *mdt, const struct mdt_library *mlib,
                  const char *filename, GError **err)
{
  int ierr;
  mod_mdt_read(&mdt->base, &mlib->base, filename, &ierr);
  if (ierr == 0) {
    mdt_setup(mdt, mlib);
    return TRUE;
  } else {
    handle_modeller_error(err);
    return FALSE;
  }
}
