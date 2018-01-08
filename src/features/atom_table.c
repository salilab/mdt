/** \file atom_table.c  Table atom feature.
 *
 *             Part of MDT, Copyright(c) 1989-2018 Andrej Sali
 */


#include "modeller.h"
#include "../mdt_index.h"
#include "../mdt_feature.h"
#include "../mdt_all_features.h"
#include "../mdt_property.h"

static int getbin(const struct mod_alignment *aln, int protein, int atom,
                  struct mdt_properties *prop,
                  const struct mdt_feature *feat,
                  const struct mdt_library *mlib,
                  const struct mod_libraries *libs, GError **err)
{
  const float *userprop = property_user(aln, protein, prop, mlib, libs,
                                        GPOINTER_TO_INT(feat->data), err);
  if (!userprop) {
    return -1;
  } else {
    return feat_to_bin(userprop[atom], feat);
  }
}

int mdt_feature_atom_table(struct mdt_library *mlib, gboolean pos2,
                           const char *table_name,
                           mdt_cb_get_property get_property,
                           gpointer data, GDestroyNotify freefunc)
{
  int iprop, ifeat;
  char *name;
  iprop = mdt_library_add_user_property(mlib, get_property, data, freefunc);

  name = g_strdup_printf("Table of atom %s", table_name);
  ifeat = mdt_feature_atom_add(mlib, name, MOD_MDTC_NONE,
                               pos2, getbin, GINT_TO_POINTER(iprop), NULL);
  g_free(name);
  return ifeat;
}
