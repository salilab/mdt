/** \file mdt_write_hdf5.c   Functions to write MDTs to binary (HDF5) files.
 *
 *             Part of MDT, Copyright(c) 1989-2007 Andrej Sali
 */

#include <glib.h>
#include "modeller.h"
#include "mdt_hdf5.h"
#include "mdt.h"
#include "util.h"

/** Write MDT data to an HDF5 file. Return TRUE on success. */
static gboolean write_mdt_data(hid_t file_id, const struct mod_mdt *mdt,
                               GError **err)
{
  herr_t ret;
  hsize_t *dims;
  int *ifeat, *istart, *iend;
  int i;

  dims = g_malloc(mdt->nfeat * sizeof(hsize_t));
  ifeat = g_malloc(mdt->nfeat * sizeof(int));
  istart = g_malloc(mdt->nfeat * sizeof(int));
  iend = g_malloc(mdt->nfeat * sizeof(int));
  for (i = 0; i < mdt->nfeat; i++) {
    const struct mod_mdt_feature *feat = &mdt->features[i];
    dims[i] = feat->nbins;
    ifeat[i] = feat->ifeat;
    istart[i] = feat->istart;
    iend[i] = feat->iend;
  }

  ret = H5LTmake_dataset_double(file_id, "/mdt", mdt->nfeat, dims, mdt->bin);
  if (ret >= 0) {
    hsize_t featdim = mdt->nfeat;
    if (H5LTmake_dataset_int(file_id, "/features", 1, &featdim, ifeat) < 0
        || H5LTmake_dataset_int(file_id, "/istart", 1, &featdim, istart) < 0
        || H5LTmake_dataset_int(file_id, "/iend", 1, &featdim, iend) < 0) {
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
  g_free(istart);
  g_free(iend);

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
