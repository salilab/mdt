/** \file mdt_write_hdf5.c   Functions to write MDTs to binary (HDF5) files.
 *
 *             Part of MDT, Copyright(c) 1989-2007 Andrej Sali
 */

#include <glib.h>
#include "modeller.h"
#include "mod_hdf5.h"
#include "mdt.h"
#include "util.h"

/** Walk the HDF5 error stack, and get the topmost error */
static herr_t errwalkfunc(int n, H5E_error_t *err_desc, void *data)
{
  GError **err = (GError **)data;
  if (n == 0) {
    g_set_error(err, MDT_ERROR, MDT_ERROR_FAILED, "%s",
                (err_desc->desc && strlen(err_desc->desc) > 0) ?
                err_desc->desc : "Generic HDF5 error");
  }
  return 0;
}

/** Convert the HDF5 error into a GError */
static void handle_hdf5_error(GError **err)
{
  H5Ewalk(H5E_WALK_DOWNWARD, errwalkfunc, err);
  if (err && *err == NULL) {
    g_set_error(err, MDT_ERROR, MDT_ERROR_FAILED, "Generic HDF5 error");
  }
}

/** Convert an HDF5 or Modeller error into a GError */
static void handle_modeller_hdf5_error(GError **err)
{
  GError *moderr = mod_error_get();
  if (moderr) {
    g_set_error(err, MDT_ERROR, MDT_ERROR_FAILED, moderr->message);
    g_error_free(moderr);
    mod_error_clear();
  } else {
    handle_hdf5_error(err);
  }
}

/** Create an HDF5 file. Return value is negative on failure. */
static hid_t mdt_hdf_create(const char *filename, unsigned flags,
                            hid_t create_plist, hid_t access_plist,
                            struct mod_file *file_info, GError **err)
{
  hid_t retval = mod_hdf_create(filename, flags, create_plist, access_plist,
                                file_info);
  if (retval < 0) {
    handle_modeller_hdf5_error(err);
  }
  return retval;
}

/** Close an HDF5 file. Return TRUE on success. */
static gboolean mdt_hdf_close(hid_t file_id, struct mod_file *file_info,
                              GError **err)
{
  if (mod_hdf_close(file_id, file_info) < 0) {
    handle_modeller_hdf5_error(err);
    return FALSE;
  } else {
    return TRUE;
  }
}

/** Write MDT data to an HDF5 file. Return TRUE on success. */
static gboolean write_mdt_data(hid_t file_id, const struct mod_mdt *mdt,
                               GError **err)
{
  herr_t ret;
  hsize_t *dims;
  int *ifeat;
  int i;

  dims = g_malloc(mdt->nfeat * sizeof(hsize_t));
  ifeat = g_malloc(mdt->nfeat * sizeof(int));
  for (i = 0; i < mdt->nfeat; i++) {
    const struct mod_mdt_feature *feat = &mdt->features[i];
    dims[i] = feat->nbins;
    ifeat[i] = feat->ifeat;
  }

  ret = H5LTmake_dataset_double(file_id, "/mdt", mdt->nfeat, dims, mdt->bin);
  if (ret >= 0) {
    hsize_t featdim[2];
    featdim[0] = mdt->nfeat;
    ret = H5LTmake_dataset_int(file_id, "/features", 1, featdim, ifeat);
  }
  g_free(dims);
  g_free(ifeat);

  if (ret < 0) {
    handle_hdf5_error(err);
    return FALSE;
  } else {
    return TRUE;
  }
}

/** Write out an MDT in HDF5 format. Return TRUE on success. */
gboolean mdt_write_hdf5(const struct mod_mdt *mdt, const char *filename,
                        GError **err)
{
  hid_t file_id;
  struct mod_file file_info;

  file_id = mdt_hdf_create(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT,
                           &file_info, err);
  if (file_id < 0 || !write_mdt_data(file_id, mdt, err)
      || !mdt_hdf_close(file_id, &file_info, err)) {
    return FALSE;
  } else {
    return TRUE;
  }
}
