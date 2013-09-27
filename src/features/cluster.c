/** \file cluster.c     Cluster feature.
 *
 *             Part of MDT, Copyright(c) 1989-2013 Andrej Sali
 */

#include "modeller.h"
#include "../mdt_index.h"
#include "../mdt_property.h"
#include "../mdt_feature.h"
#include "../mdt_all_features.h"

struct feature_data {
  /* Keys are (bin1,bin2) pairs (of the child features); values are the
     corresponding bin index of the cluster feature */
  GHashTable *map;
};

static int getbin(int bin1, int bin2, struct mdt_properties *prop,
                  const struct mdt_feature *feat,
                  const struct mdt_library *mlib,
                  const struct mod_libraries *libs, GError **err)
{
  struct feature_data *feat_data = (struct feature_data *)feat->data;
  gpointer *val = g_hash_table_lookup(feat_data->map,
                                      MAKE_HASH_KEY_ASYMMETRIC(bin1, bin2));
  if (val) {
    return GPOINTER_TO_INT(val);
  } else {
    return mdt_feature_undefined_bin_get(feat);
  }
}

static void free_data(void *data)
{
  struct feature_data *feat_data = (struct feature_data *)data;
  g_hash_table_destroy(feat_data->map);
}

void mdt_cluster_add(struct mdt_feature *feat, int bin1, int bin2, int bin)
{
  struct feature_data *feat_data = (struct feature_data *)feat->data;
  g_hash_table_insert(feat_data->map,
                      MAKE_HASH_KEY_ASYMMETRIC(bin1, bin2),
                      GINT_TO_POINTER(bin));
}

int mdt_feature_cluster(struct mdt_library *mlib, int ifeat1, int ifeat2,
                        GError **err)
{
  struct feature_data *feat_data;
  int ifeat;

  feat_data = g_malloc(sizeof(struct feature_data));
  feat_data->map = g_hash_table_new(NULL, NULL);
  ifeat = mdt_feature_group_add(mlib, "Cluster", MOD_MDTC_NONE, ifeat1, ifeat2,
                                getbin, feat_data, free_data, err);
  if (ifeat < 0) {
    free_data(feat_data);
  }
  return ifeat;
}
