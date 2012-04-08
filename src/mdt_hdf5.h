/** \file mdt_hdf5.h      Functions to handle HDF5 files.
 *
 *             Part of MDT, Copyright(c) 1989-2012 Andrej Sali
 */

#ifndef __MDT_HDF5_H
#define __MDT_HDF5_H

#include <glib.h>
#include "mdt_config.h"
#include "mod_hdf5.h"

G_BEGIN_DECLS

/** Convert the HDF5 error into a GError */
MDTDLLLOCAL
void handle_hdf5_error(GError **err);

/** Convert an HDF5 or Modeller error into a GError */
MDTDLLLOCAL
void handle_modeller_hdf5_error(GError **err);

/** Create an HDF5 file. Return value is negative on failure. */
MDTDLLLOCAL
hid_t mdt_hdf_create(const char *filename, unsigned flags,
                     hid_t create_plist, hid_t access_plist,
                     struct mod_file *file_info, GError **err);

/** Open an existing HDF5 file. Return value is negative on failure. */
MDTDLLLOCAL
hid_t mdt_hdf_open(const char *filename, unsigned flags, hid_t access_plist,
                   struct mod_file *file_info, GError **err);

/** Close an HDF5 file. Return TRUE on success. */
MDTDLLLOCAL
gboolean mdt_hdf_close(hid_t file_id, struct mod_file *file_info, GError **err);

#endif  /*  __MDT_HDF5_H */
