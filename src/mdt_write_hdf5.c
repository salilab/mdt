/** \file mdt_write_hdf5.c   Functions to write MDTs to binary (HDF5) files.
 *
 *             Part of MDT, Copyright(c) 1989-2013 Andrej Sali
 */

#include <glib.h>
#include "modeller.h"
#include "mdt_hdf5.h"
#include "mdt_feature.h"
#include "mdt.h"
#include "util.h"

/** Write a float attribute. Return TRUE on success. */
static gboolean write_float_attribute(hid_t loc_id, const char *name,
                                      float value)
{
  hid_t attr;
  return (attr = H5Acreate(loc_id, name, H5T_NATIVE_FLOAT, H5S_SCALAR,
                           H5P_DEFAULT, H5P_DEFAULT)) >= 0
         && H5Awrite(attr, H5T_NATIVE_FLOAT, &value) >= 0
         && H5Aclose(attr) >= 0;
}

static gboolean write_ifeat(hid_t group_id, const struct mdt *mdt,
                            const struct mdt_library *mlib, int ifeat,
                            GHashTable *seen_lib_funcs)
{
  char *group_name;
  hid_t featgroup_id;
  const struct mod_mdt_libfeature *feat = &mlib->base.features[ifeat - 1];
  const struct mdt_feature *mfeat = &g_array_index(mlib->features,
                                                   struct mdt_feature,
                                                   ifeat - 1);

  group_name = g_strdup_printf("feature%d", ifeat);
  featgroup_id = H5Gcreate(group_id, group_name, H5P_DEFAULT, H5P_DEFAULT,
                           H5P_DEFAULT);
  g_free(group_name);
  if (featgroup_id < 0) {
    return FALSE;
  }

  if (mfeat->writelibfunc
      && !g_hash_table_lookup(seen_lib_funcs, mfeat->writelibfunc)) {
    /* Only call each function once per file */
    g_hash_table_insert(seen_lib_funcs, mfeat->writelibfunc,
                        GINT_TO_POINTER(1));
    /* Note that datasets go in the top-level group, not the per-feature
       group */
    if (!mfeat->writelibfunc(group_id, mfeat, mlib)) {
      return FALSE;
    }
  }

  if (mfeat->writefunc) {
    if (!mfeat->writefunc(featgroup_id, mfeat, mlib)) {
      return FALSE;
    }
  }

  if (mfeat->uniform_bins) {
    if (!write_float_attribute(featgroup_id, "bin_width", mfeat->bin_width)
        || !write_float_attribute(featgroup_id, "first_bin",
                                  feat->bins[0].rang1)) {
      return FALSE;
    }
  }

  return H5Gclose(featgroup_id) >= 0;
}

/** Write a single MDT library feature to file. Return TRUE on success. */
static gboolean write_library_feature(hid_t group_id, const struct mdt *mdt,
                                      const struct mdt_library *mlib, int nfeat,
                                      GHashTable *seen_lib_funcs)
{
  const struct mod_mdt_feature *feat = &mdt->base.features[nfeat];
  const struct mdt_feature *mfeat = &g_array_index(mlib->features,
                                                   struct mdt_feature,
                                                   feat->ifeat - 1);
  if (mfeat->type == MDT_FEATURE_GROUP) {
    return write_ifeat(group_id, mdt, mlib, mfeat->u.group.ifeat1,
                       seen_lib_funcs)
           && write_ifeat(group_id, mdt, mlib, mfeat->u.group.ifeat2,
                          seen_lib_funcs)
           && write_ifeat(group_id, mdt, mlib, feat->ifeat, seen_lib_funcs);
  } else {
    return write_ifeat(group_id, mdt, mlib, feat->ifeat, seen_lib_funcs);
  }
}

/** Write MDT library information to file. Return TRUE on success. */
static gboolean write_library_info(hid_t file_id, const struct mdt *mdt,
                                   const struct mdt_library *mlib)
{
  hid_t group_id;
  GHashTable *seen_lib_funcs;
  gboolean retval = TRUE;
  int i;

  group_id = H5Gcreate(file_id, "/library", H5P_DEFAULT, H5P_DEFAULT,
                       H5P_DEFAULT);
  if (group_id < 0) {
    return FALSE;
  }

  seen_lib_funcs = g_hash_table_new(NULL, NULL);
  for (i = 0; i < mdt->base.nfeat && retval; i++) {
    retval = retval && write_library_feature(group_id, mdt, mlib, i,
                                             seen_lib_funcs);
  }
  g_hash_table_destroy(seen_lib_funcs);

  return retval && H5Gclose(group_id) >= 0;
}

/** Write MDT feature names to file. Return TRUE on success. */
static gboolean write_mdt_feature_names(hid_t file_id, char **name,
                                        hsize_t featdim[])
{
  hid_t dtype;

  dtype = H5Tcopy(H5T_C_S1);
  if (dtype < 0) {
    return FALSE;
  } else {
    gboolean ok = TRUE;
    if (H5Tset_size(dtype, H5T_VARIABLE) < 0
        || H5Tset_strpad(dtype, H5T_STR_NULLTERM) < 0
        || H5LTmake_dataset(file_id, "/feature_names", 1, featdim, dtype,
                            name) < 0) {
      ok = FALSE;
    }
    return (H5Tclose(dtype) >= 0 && ok);
  }
}

/** Write MDT data to an HDF5 file. Return TRUE on success. */
static gboolean write_mdt_data(hid_t file_id, const struct mdt *mdt,
                               const struct mdt_library *mlib, GError **err)
{
  herr_t ret;
  hid_t type_id;
  hsize_t *dims;
  int *ifeat, *offset, *nbins;
  char **name;
  int i;

  dims = g_malloc(mdt->base.nfeat * sizeof(hsize_t));
  ifeat = g_malloc(mdt->base.nfeat * sizeof(int));
  offset = g_malloc(mdt->base.nfeat * sizeof(int));
  nbins = g_malloc(mdt->base.nfeat * sizeof(int));
  name = g_malloc(mdt->base.nfeat * sizeof(char *));
  for (i = 0; i < mdt->base.nfeat; i++) {
    const struct mod_mdt_libfeature *libfeat;
    const struct mod_mdt_feature *feat = &mdt->base.features[i];
    libfeat = &mlib->base.features[feat->ifeat - 1];
    dims[i] = feat->nbins;
    ifeat[i] = feat->ifeat;
    offset[i] = feat->istart - 1;
    nbins[i] = libfeat->nbins;
    name[i] = libfeat->name;
  }

  type_id = mdt_get_hdf5_type(&mdt->base);
  ret = H5LTmake_dataset(file_id, "/mdt", mdt->base.nfeat, dims, type_id,
                         mdt->base.bindata);
  if (ret >= 0) {
    hsize_t featdim = mdt->base.nfeat;
    if (H5LTmake_dataset_int(file_id, "/features", 1, &featdim, ifeat) < 0
        || H5LTmake_dataset_int(file_id, "/offset", 1, &featdim, offset) < 0
        || H5LTmake_dataset_int(file_id, "/nbins", 1, &featdim, nbins) < 0
        || !write_mdt_feature_names(file_id, name, &featdim)
        || !write_library_info(file_id, mdt, mlib)) {
      ret = -1;
    }
  }
  if (ret >= 0) {
    hsize_t featdim = 1;
    char is_pdf = (mdt->pdf ? 1 : 0);
    if (H5LTmake_dataset_int(file_id, "/n_alignments", 1, &featdim,
                             &mdt->nalns) < 0
        || H5LTmake_dataset_int(file_id, "/n_proteins", 1, &featdim,
                                &mdt->n_proteins) < 0
        || H5LTmake_dataset_int(file_id, "/n_protein_pairs", 1, &featdim,
                                &mdt->n_protein_pairs) < 0
        || H5LTmake_dataset_double(file_id, "/sample_size", 1, &featdim,
                                   &mdt->sample_size) < 0
        || H5LTmake_dataset_char(file_id, "/pdf", 1, &featdim, &is_pdf) < 0) {
      ret = -1;
    }
  }
  g_free(dims);
  g_free(ifeat);
  g_free(offset);
  g_free(name);
  g_free(nbins);

  if (ret < 0) {
    handle_hdf5_error(err);
    return FALSE;
  } else {
    return TRUE;
  }
}

/** Write out an MDT in HDF5 format. Return TRUE on success. */
gboolean mdt_write_hdf5(const struct mdt *mdt, const struct mdt_library *mlib,
                        const char *filename, GError **err)
{
  hid_t file_id;
  struct mod_file file_info;

  file_id = mdt_hdf_create(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT,
                           &file_info, err);
  if (file_id < 0 || !write_mdt_data(file_id, mdt, mlib, err)
      || !mdt_hdf_close(file_id, &file_info, err)) {
    return FALSE;
  } else {
    return TRUE;
  }
}
