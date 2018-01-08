/** \file cluster.c     Cluster feature.
 *
 *             Part of MDT, Copyright(c) 1989-2018 Andrej Sali
 */

#include "modeller.h"
#include "../mdt_error.h"
#include "../mdt_index.h"
#include "../mdt_property.h"
#include "../mdt_feature.h"
#include "../mdt_hdf5.h"
#include "../mdt_all_features.h"

struct feature_data {
  /* Keys are (bin1,bin2) pairs (of the child features); values are the
     corresponding bin index of the cluster feature */
  GHashTable *map;
};

struct store_cluster_data {
  int *cluster;
};

void store_cluster(gpointer key, gpointer value, gpointer user_data)
{
  struct store_cluster_data *data = (struct store_cluster_data *)user_data;
  *(data->cluster++) = GET_HASH_KEY_HIGH(key);
  *(data->cluster++) = GET_HASH_KEY_LOW(key);
  *(data->cluster++) = GPOINTER_TO_INT(value);
}

static gboolean writefunc(hid_t loc_id, const struct mdt_feature *feat,
                          const struct mdt_library *mlib)
{
  struct store_cluster_data user_data;
  int *cluster;
  herr_t retval;
  struct feature_data *feat_data = (struct feature_data *)feat->data;
  hsize_t dims[2];
  dims[0] = g_hash_table_size(feat_data->map);
  dims[1] = 3;

  cluster = g_malloc(dims[0] * dims[1] * sizeof(int));
  user_data.cluster = cluster;
  g_hash_table_foreach(feat_data->map, store_cluster, &user_data);
  retval = H5LTmake_dataset_int(loc_id, "clusters", 2, dims, cluster);
  g_free(cluster);
  return retval >= 0;
}

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

static gboolean check_bin(const struct mod_mdt_libfeature *feat, int bin,
                          int nfeat, GError **err)
{
  if (bin < 0 || bin >= feat->nbins) {
    g_set_error(err, MDT_ERROR, MDT_ERROR_VALUE,
                "Bin %d index (%d) is out of range 0-%d",
                nfeat, bin, feat->nbins - 1);
    return FALSE;
  } else {
    return TRUE;
  }
}

gboolean mdt_cluster_add(struct mdt_library *mlib, int ifeat,
                         int bin1, int bin2, int bin, GError **err)
{
  struct mdt_feature *mfeat;
  struct mod_mdt_libfeature *feat, *feat1, *feat2;
  struct feature_data *feat_data;

  /* Note that we take bin indexes starting from 0 (for consistency with
     the rest of the Python interface) but internally they are stored
     starting from 1 */

  feat = &mlib->base.features[ifeat - 1];
  mfeat = &g_array_index(mlib->features, struct mdt_feature, ifeat - 1);
  if (mfeat->type != MDT_FEATURE_GROUP || mfeat->freefunc != free_data) {
    g_set_error(err, MDT_ERROR, MDT_ERROR_VALUE,
                "Feature is not a cluster feature");
    return FALSE;
  }
  feat1 = &mlib->base.features[mfeat->u.group.ifeat1 - 1];
  feat2 = &mlib->base.features[mfeat->u.group.ifeat2 - 1];
  if (!check_bin(feat1, bin1, 1, err) || !check_bin(feat2, bin2, 2, err)) {
    return FALSE;
  }

  if (bin < 0 || bin >= feat->nbins) {
    g_set_error(err, MDT_ERROR, MDT_ERROR_VALUE,
                "Output bin index (%d) is out of range 0-%d",
                bin, feat->nbins - 1);
    return FALSE;
  }

  feat_data = (struct feature_data *)mfeat->data;
  g_hash_table_insert(feat_data->map,
                      MAKE_HASH_KEY_ASYMMETRIC(bin1 + 1, bin2 + 1),
                      GINT_TO_POINTER(bin + 1));
  return TRUE;
}

int mdt_feature_cluster(struct mdt_library *mlib, int ifeat1, int ifeat2,
                        int nbins, GError **err)
{
  struct feature_data *feat_data;
  int ifeat, i;

  feat_data = g_malloc(sizeof(struct feature_data));
  feat_data->map = g_hash_table_new(NULL, NULL);
  ifeat = mdt_feature_group_add(mlib, "Cluster", MOD_MDTC_NONE, ifeat1, ifeat2,
                                getbin, feat_data, free_data, err);
  if (ifeat < 0) {
    free_data(feat_data);
  } else {
    struct mod_mdt_libfeature *feat = &mlib->base.features[ifeat - 1];
    mdt_feature_set_write_callback(mlib, ifeat, writefunc);
    mod_mdt_libfeature_nbins_set(feat, nbins + 1);
    for (i = 0; i < nbins; ++i) {
      g_free(feat->bins[i].symbol);
      feat->bins[i].symbol = g_strdup("C");
      feat->bins[i].rang1 = i;
      feat->bins[i].rang2 = i + 1;
    }
    g_free(feat->bins[nbins].symbol);
    feat->bins[nbins].symbol = g_strdup("U");
  }
  return ifeat;
}
