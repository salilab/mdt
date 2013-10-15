/** \file mdt_write_hdf5.c   Functions to write MDTs to binary (HDF5) files.
 *
 *             Part of MDT, Copyright(c) 1989-2013 Andrej Sali
 */

#include <glib.h>
#include <math.h>
#include "modeller.h"
#include "mdt_hdf5.h"
#include "mdt_feature.h"
#include "mdt.h"
#include "util.h"

static gboolean write_ifeat(hid_t group_id, const struct mdt *mdt,
                            const struct mdt_library *mlib, int ifeat)
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

  if (mfeat->writefunc) {
    if (!mfeat->writefunc(featgroup_id, mfeat, mlib)) {
      return FALSE;
    }
  }

  if (mfeat->uniform_bins) {
    if (!mdt_hdf5_write_float_attr(featgroup_id, "bin_width", 1,
                                   &mfeat->bin_width)
        || !mdt_hdf5_write_float_attr(featgroup_id, "first_bin", 1,
                                      &feat->bins[0].rang1)) {
      return FALSE;
    }
  }

  return H5Gclose(featgroup_id) >= 0;
}

struct callback_data {
  gboolean retval;
  hid_t loc_id;
  const struct mdt_library *mlib;
};

static void write_callbacks(gpointer key, gpointer value, gpointer user_data)
{
  mdt_cb_write_lib writelibfunc = (mdt_cb_write_lib)key;
  struct callback_data *data = (struct callback_data *)user_data;

  if (data->retval) {
    data->retval = writelibfunc(data->loc_id, data->mlib);
  }
}

/** Write all libraries used by this MDT. Return TRUE on success. */
static gboolean write_used_libs(hid_t group_id, const struct mdt *mdt,
                                const struct mdt_library *mlib)
{
  struct callback_data data;
  data.retval = TRUE;
  data.loc_id = group_id;
  data.mlib = mlib;
  g_hash_table_foreach(mdt->write_lib_funcs, write_callbacks, &data);
  return data.retval;
}

/** Write a single MDT library feature to file. Return TRUE on success. */
static gboolean write_library_feature(hid_t group_id, const struct mdt *mdt,
                                      const struct mdt_library *mlib, int nfeat)
{
  const struct mod_mdt_feature *feat = &mdt->base.features[nfeat];
  const struct mdt_feature *mfeat = &g_array_index(mlib->features,
                                                   struct mdt_feature,
                                                   feat->ifeat - 1);
  if (mfeat->type == MDT_FEATURE_GROUP) {
    return write_ifeat(group_id, mdt, mlib, mfeat->u.group.ifeat1)
           && write_ifeat(group_id, mdt, mlib, mfeat->u.group.ifeat2)
           && write_ifeat(group_id, mdt, mlib, feat->ifeat);
  } else {
    return write_ifeat(group_id, mdt, mlib, feat->ifeat);
  }
}

/** Write MDT scan information to file. Return TRUE on success. */
static gboolean write_scan_info(hid_t file_id,
                                const struct mdt_scan_parameters *params)
{
  hid_t group_id;
  /* Don't do anything if no scan information is available */
  if (!params->scan_called) {
    return TRUE;
  }
  return (group_id = H5Gcreate(file_id, "/scan", H5P_DEFAULT, H5P_DEFAULT,
                               H5P_DEFAULT)) >= 0
         && mdt_hdf5_write_int_attr(group_id, "residue_span_range", 4,
                                    params->residue_span_range)
         && mdt_hdf5_write_int_attr(group_id, "chain_span_range", 4,
                                    params->chain_span_range)
         && mdt_hdf5_write_int_attr(group_id, "bond_span_range", 2,
                                    params->bond_span_range)
         && mdt_hdf5_write_int_attr(group_id, "disulfide", 1,
                                    &params->disulfide)
         && mdt_hdf5_write_int_attr(group_id, "exclude_bonds", 1,
                                    &params->exclude_bonds)
         && mdt_hdf5_write_int_attr(group_id, "exclude_angles", 1,
                                    &params->exclude_angles)
         && mdt_hdf5_write_int_attr(group_id, "exclude_dihedrals", 1,
                                    &params->exclude_dihedrals)
         && mdt_hdf5_write_int_attr(group_id, "sympairs", 1,
                                    &params->sympairs)
         && mdt_hdf5_write_int_attr(group_id, "symtriples", 1,
                                    &params->symtriples)
         && mdt_hdf5_write_float_attr(group_id, "distngh", 1,
                                      &params->distngh)
         && mdt_hdf5_write_int_attr(group_id, "surftyp", 1,
                                    &params->surftyp)
         && mdt_hdf5_write_int_attr(group_id, "accessibility_type", 1,
                                    &params->accessibility_type)
         && H5Gclose(group_id) >= 0;
}

/** Write MDT library information to file. Return TRUE on success. */
static gboolean write_library_info(hid_t file_id, const struct mdt *mdt,
                                   const struct mdt_library *mlib)
{
  hid_t group_id;
  gboolean retval = TRUE;
  int i;

  group_id = H5Gcreate(file_id, "/library", H5P_DEFAULT, H5P_DEFAULT,
                       H5P_DEFAULT);
  if (group_id < 0) {
    return FALSE;
  }

  retval = retval && write_used_libs(group_id, mdt, mlib);
  for (i = 0; i < mdt->base.nfeat && retval; i++) {
    retval = retval && write_library_feature(group_id, mdt, mlib, i);
  }

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

/** Determine the dimensions of a single chunk. */
static void get_chunk_dims(const hsize_t *dims, hsize_t *chunk_dims,
                           int ndims, int chunk_size)
{
  int i;
  hsize_t total_size;
  double div;

  /* First guess: assume dimensions are similar */
  total_size = 1;
  for (i = 0; i < ndims; i++) {
    total_size *= dims[i];
  }
  div = pow((double)total_size / (double)chunk_size, 1. / (double)ndims);

  for (i = 0; i < ndims; i++) {
    chunk_dims[i] = (int)(dims[i] / div);
    chunk_dims[i] = CLAMP(chunk_dims[i], 1, dims[i]);
  }
}

/** Write the MDT table itself to file. Return TRUE on success. */
static gboolean write_table(hid_t file_id, const struct mdt *mdt, int gzip,
                            int chunk_size)
{
  gboolean ret;
  hsize_t *dims, *chunk_dims = NULL;
  hid_t type_id, prop_id, space_id, set_id;
  int i;
  dims = g_malloc(mdt->base.nfeat * sizeof(hsize_t));
  for (i = 0; i < mdt->base.nfeat; i++) {
    const struct mod_mdt_feature *feat = &mdt->base.features[i];
    dims[i] = feat->nbins;
  }
  type_id = mdt_get_hdf5_type(&mdt->base);

  ret = TRUE;
  if (gzip >= 0) {
    chunk_dims = g_malloc(mdt->base.nfeat * sizeof(hsize_t));
    get_chunk_dims(dims, chunk_dims, mdt->base.nfeat, chunk_size);
    ret = (prop_id = H5Pcreate(H5P_DATASET_CREATE)) >= 0
          && H5Pset_chunk(prop_id, mdt->base.nfeat, chunk_dims) >= 0
          && H5Pset_deflate(prop_id, gzip) >= 0;
  } else {
    prop_id = H5P_DEFAULT;
  }

  ret = ret
        && (space_id = H5Screate_simple(mdt->base.nfeat, dims, NULL)) >= 0
        && (set_id = H5Dcreate(file_id, "/mdt", type_id, space_id, H5P_DEFAULT,
                               prop_id, H5P_DEFAULT)) >= 0
        && H5Dwrite(set_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                    mdt->base.bindata) >= 0
        && H5Dclose(set_id) >= 0
        && H5Sclose(space_id) >= 0;

  g_free(dims);
  if (gzip >= 0) {
    g_free(chunk_dims);
    ret = ret && H5Pclose(prop_id) >= 0;
  }
  return ret;
}

/** Write MDT data to an HDF5 file. Return TRUE on success. */
static gboolean write_mdt_data(hid_t file_id, const struct mdt *mdt,
                               const struct mdt_library *mlib, int gzip,
                               int chunk_size, GError **err)
{
  herr_t ret = 0;
  hsize_t featdim;
  int *ifeat, *offset, *nbins;
  char **name;
  int i;

  ifeat = g_malloc(mdt->base.nfeat * sizeof(int));
  offset = g_malloc(mdt->base.nfeat * sizeof(int));
  nbins = g_malloc(mdt->base.nfeat * sizeof(int));
  name = g_malloc(mdt->base.nfeat * sizeof(char *));
  for (i = 0; i < mdt->base.nfeat; i++) {
    const struct mod_mdt_libfeature *libfeat;
    const struct mod_mdt_feature *feat = &mdt->base.features[i];
    libfeat = &mlib->base.features[feat->ifeat - 1];
    ifeat[i] = feat->ifeat;
    offset[i] = feat->istart - 1;
    nbins[i] = libfeat->nbins;
    name[i] = libfeat->name;
  }

  featdim = mdt->base.nfeat;
  if (!write_table(file_id, mdt, gzip, chunk_size)
      || H5LTmake_dataset_int(file_id, "/features", 1, &featdim, ifeat) < 0
      || H5LTmake_dataset_int(file_id, "/offset", 1, &featdim, offset) < 0
      || H5LTmake_dataset_int(file_id, "/nbins", 1, &featdim, nbins) < 0
      || !write_mdt_feature_names(file_id, name, &featdim)
      || !write_scan_info(file_id, &mdt->scan_params)
      || !write_library_info(file_id, mdt, mlib)) {
    ret = -1;
  }

  if (ret >= 0) {
    char is_pdf = (mdt->pdf ? 1 : 0);
    featdim = 1;
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
                        const char *filename, int gzip, int chunk_size,
                        GError **err)
{
  hid_t file_id;
  struct mod_file file_info;

  file_id = mdt_hdf_create(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT,
                           &file_info, err);
  if (file_id < 0 || !write_mdt_data(file_id, mdt, mlib, gzip, chunk_size, err)
      || !mdt_hdf_close(file_id, &file_info, err)) {
    return FALSE;
  } else {
    return TRUE;
  }
}
