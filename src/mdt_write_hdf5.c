/** \file mdt_write_hdf5.c   Functions to write MDTs to binary (HDF5) files.
 *
 *             Part of MDT, Copyright(c) 1989-2013 Andrej Sali
 */

#include <glib.h>
#include "modeller.h"
#include "mdt_hdf5.h"
#include "mdt.h"
#include "util.h"

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
        || write_mdt_feature_names(file_id, name, &featdim) < 0) {
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
