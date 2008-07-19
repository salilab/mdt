/** \file mdt_add.c       Functions to add an MDT to another.
 *
 *             Part of MDT, Copyright(c) 1989-2008 Andrej Sali
 */

#include "mdt.h"

/** Check to see if two MDTs are compatible (same features, size, etc.).
    \return NULL, or a string describing the incompatibility. */
static char *check_mdts_compatible(const struct mdt *mdt1,
                                   const struct mdt *mdt2)
{
  if (mdt1->base.nfeat != mdt2->base.nfeat) {
    return g_strdup_printf("incompatible number of features (%d vs. %d)",
                           mdt1->base.nfeat, mdt2->base.nfeat);
  } else {
    int i;
    for (i = 0; i < mdt1->base.nfeat; ++i) {
      struct mod_mdt_feature *feat1 = &mdt1->base.features[i];
      struct mod_mdt_feature *feat2 = &mdt2->base.features[i];
      if (feat1->ifeat != feat2->ifeat) {
        return g_strdup_printf("feature %d types are incompatible (%d vs. %d)",
                               i, feat1->ifeat, feat2->ifeat);
      } else if (feat1->istart != feat2->istart) {
        return g_strdup_printf("feature %d starts are incompatible (%d vs. %d)",
                               i, feat1->istart, feat2->istart);
      } else if (feat1->nbins != feat2->nbins) {
        return g_strdup_printf("feature %d nbins are incompatible (%d vs. %d)",
                               i, feat1->nbins, feat2->nbins);
      }
    }
  }
  return NULL;
}

/** Add mdt2 into mdt1. Return TRUE on success. */
gboolean mdt_add(struct mdt *mdt1, const struct mdt *mdt2, GError **err)
{
  size_t i;
  char *msg = check_mdts_compatible(mdt1, mdt2);
  if (msg) {
    g_set_error(err, MDT_ERROR, MDT_ERROR_VALUE, "Cannot add MDTs: %s", msg);
    g_free(msg);
    return FALSE;
  }
  for (i = 0; i < mdt1->base.nelems; ++i) {
    mod_mdt_bin_set(&mdt1->base, i, mod_mdt_bin_get(&mdt1->base, i)
                                    + mod_mdt_bin_get(&mdt2->base, i));
  }
  return TRUE;
}
