/** \file mdt_hdf5.c      Functions to handle HDF5 files.
 *
 *             Part of MDT, Copyright(c) 1989-2021 Andrej Sali
 */

#include <glib.h>
#include <string.h>
#include "mod_error.h"
#include "mdt_hdf5.h"
#include "mdt_error.h"

/** Walk the HDF5 error stack, and get the topmost error */
static herr_t errwalkfunc(unsigned n, const H5E_error2_t *err_desc, void *data)
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
void handle_hdf5_error(GError **err)
{
  H5Ewalk2(H5E_DEFAULT, H5E_WALK_DOWNWARD, errwalkfunc, err);
  if (err && *err == NULL) {
    g_set_error(err, MDT_ERROR, MDT_ERROR_FAILED, "Generic HDF5 error");
  }
}

/** Convert an HDF5 or Modeller error into a GError */
void handle_modeller_hdf5_error(GError **err)
{
  GError *moderr = mod_error_get();
  if (moderr) {
    g_set_error(err, MDT_ERROR, MDT_ERROR_FAILED, "%s", moderr->message);
    g_error_free(moderr);
    mod_error_clear();
  } else {
    handle_hdf5_error(err);
  }
}

/** Create an HDF5 file. Return value is negative on failure. */
hid_t mdt_hdf_create(const char *filename, unsigned flags,
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

/** Open an existing HDF5 file. Return value is negative on failure. */
hid_t mdt_hdf_open(const char *filename, unsigned flags, hid_t access_plist,
                   struct mod_file *file_info, GError **err)
{
  hid_t retval = mod_hdf_open(filename, flags, access_plist, file_info);
  if (retval < 0) {
    handle_modeller_hdf5_error(err);
  }
  return retval;
}

/** Close an HDF5 file. Return TRUE on success. */
gboolean mdt_hdf_close(hid_t file_id, struct mod_file *file_info, GError **err)
{
  if (mod_hdf_close(file_id, file_info) < 0) {
    handle_modeller_hdf5_error(err);
    return FALSE;
  } else {
    return TRUE;
  }
}

/** Write a float attribute. Return TRUE on success. */
gboolean mdt_hdf5_write_float_attr(hid_t loc_id, const char *name,
                                   int dim, const float *value)
{
  hid_t attr, dataspace_id;
  hsize_t dims = dim;
  return (dataspace_id = H5Screate_simple(1, &dims, NULL)) >= 0
         && (attr = H5Acreate2(loc_id, name, H5T_NATIVE_FLOAT, dataspace_id,
                               H5P_DEFAULT, H5P_DEFAULT)) >= 0
         && H5Awrite(attr, H5T_NATIVE_FLOAT, value) >= 0
         && H5Aclose(attr) >= 0
         && H5Sclose(dataspace_id) >= 0;
}

/** Write an int attribute. Return TRUE on success. */
gboolean mdt_hdf5_write_int_attr(hid_t loc_id, const char *name,
                                 int dim, const int *value)
{
  hid_t attr, dataspace_id;
  hsize_t dims = dim;
  return (dataspace_id = H5Screate_simple(1, &dims, NULL)) >= 0
         && (attr = H5Acreate2(loc_id, name, H5T_NATIVE_INT, dataspace_id,
                               H5P_DEFAULT, H5P_DEFAULT)) >= 0
         && H5Awrite(attr, H5T_NATIVE_INT, value) >= 0
         && H5Aclose(attr) >= 0
         && H5Sclose(dataspace_id) >= 0;
}
