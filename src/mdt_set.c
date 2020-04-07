/** \file mdt_set.c        Functions to set elements in MDTs.
 *
 *             Part of MDT, Copyright(c) 1989-2020 Andrej Sali
 */

#include "mdt.h"
#include "util.h"

/** Set an element in an MDT. Return TRUE on success. */
gboolean mdt_set(struct mod_mdt *mdt, const int indices[], int n_indices,
                 double val, GError **err)
{
  int bin_index;
  if (get_bin_index(mdt, indices, n_indices, &bin_index, err)) {
    mod_mdt_bin_set(mdt, bin_index, val);
    return TRUE;
  } else {
    return FALSE;
  }
}
