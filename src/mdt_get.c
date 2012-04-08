/** \file mdt_get.c        Functions to get elements from MDTs.
 *
 *             Part of MDT, Copyright(c) 1989-2012 Andrej Sali
 */

#include "mdt.h"
#include "util.h"

/** Get an element from an MDT. */
double mdt_get(const struct mod_mdt *mdt, const int indices[], int n_indices,
               GError **err)
{
  int bin_index;
  if (get_bin_index(mdt, indices, n_indices, &bin_index, err)) {
    return mod_mdt_bin_get(mdt, bin_index);
  } else {
    return 0.0;
  }
}
