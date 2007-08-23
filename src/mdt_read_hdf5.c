/** \file mdt_read_hdf5.c   Functions to read MDTs from binary (HDF5) files.
 *
 *             Part of MDT, Copyright(c) 1989-2007 Andrej Sali
 */

#include <glib.h>
#include "modeller.h"
#include "mdt_hdf5.h"
#include "mdt.h"
#include "util.h"

/** Read MDT feature information from an HDF5 file. Return TRUE on success. */
static gboolean read_mdt_features(hid_t file_id, struct mod_mdt *mdt,
                                  GError **err)
{
  hsize_t nfeat;
  int *ifeat, *istart, *iend;
  herr_t ret;

  ret = mod_dataset_get_ndsize(file_id, "/features", 1, &nfeat);
  if (ret < 0) {
    handle_modeller_hdf5_error(err);
    return FALSE;
  }
  mod_mdt_nfeat_set(mdt, nfeat);
  ifeat = g_malloc(mdt->nfeat * sizeof(int));
  istart = g_malloc(mdt->nfeat * sizeof(int));
  iend = g_malloc(mdt->nfeat * sizeof(int));
  if (mod_dataset_read_int(file_id, "/features", 1, &nfeat, ifeat) >= 0
      && mod_dataset_read_int(file_id, "/istart", 1, &nfeat, istart) >= 0
      && mod_dataset_read_int(file_id, "/iend", 1, &nfeat, iend) >= 0) {
    int i;
    for (i = 0; i < mdt->nfeat; i++) {
      struct mod_mdt_feature *feat = &mdt->features[i];
      feat->ifeat = ifeat[i];
      feat->istart = istart[i];
      feat->iend = iend[i];
      feat->nbins = iend[i] - istart[i] + 1;
    }
  } else {
    handle_modeller_hdf5_error(err);
    ret = -1;
  }
  g_free(ifeat);
  g_free(istart);
  g_free(iend);
  return (ret >= 0);
}


/** Check MDT feature information for sanity. */
static gboolean check_mdt_features(const struct mod_mdt *mdt,
                                   const struct mdt_library *mlib,
                                   GError **err)
{
  int i;
  for (i = 0; i < mdt->nfeat; i++) {
    const struct mod_mdt_libfeature *libfeat;
    const struct mod_mdt_feature *feat = &mdt->features[i];

    if (!check_feature_type(feat->ifeat, mlib, err)) {
      return FALSE;
    }

    libfeat = &mlib->base.features[feat->ifeat - 1];
    if (feat->nbins > libfeat->nbins) {
      g_set_error(err, MDT_ERROR, MDT_ERROR_FAILED,
                  "Number of bins in MDT table (%d) is greater than "
                  "that in bin file (%d) for feature number %d, feature "
                  "type %d", feat->nbins, libfeat->nbins, i, feat->ifeat);
      return FALSE;
    }
  }
  return TRUE;
}


/** Read MDT data from an HDF5 file. Return TRUE on success. */
static gboolean read_mdt_data(hid_t file_id, struct mod_mdt *mdt, GError **err)
{
  herr_t ret;
  hsize_t *dims;
  int i, nelems;

  dims = g_malloc(mdt->nfeat * sizeof(hsize_t));
  nelems = 1;
  for (i = 0; i < mdt->nfeat; i++) {
    const struct mod_mdt_feature *feat = &mdt->features[i];
    dims[i] = feat->nbins;
    nelems *= feat->nbins;
  }

  mod_mdt_nelems_set(mdt, nelems);
  ret = mod_dataset_read_double(file_id, "/mdt", mdt->nfeat, dims, mdt->bin);
  if (ret >= 0) {
    hsize_t featdim = 1;
    char is_pdf;
    if (mod_dataset_read_int(file_id, "/n_alignments", 1, &featdim,
                             &mdt->nalns) < 0
        || mod_dataset_read_int(file_id, "/n_proteins", 1, &featdim,
                                &mdt->n_proteins) < 0
        || mod_dataset_read_int(file_id, "/n_protein_pairs", 1, &featdim,
                                &mdt->n_protein_pairs) < 0
        || mod_dataset_read_double(file_id, "/sample_size", 1, &featdim,
                                   &mdt->sample_size) < 0
        || mod_dataset_read_char(file_id, "/pdf", 1, &featdim, &is_pdf) < 0) {
      ret = -1;
    } else {
      mdt->pdf = (is_pdf == 1);
    }
  }
  g_free(dims);

  if (ret < 0) {
    handle_modeller_hdf5_error(err);
    return FALSE;
  } else {
    return TRUE;
  }
}

/** Read in an MDT in HDF5 format. Return TRUE on success. */
gboolean mdt_read_hdf5(struct mod_mdt *mdt, const struct mdt_library *mlib,
                       const char *filename, GError **err)
{
  hid_t file_id;
  struct mod_file file_info;

  file_id = mdt_hdf_open(filename, H5F_ACC_RDONLY, H5P_DEFAULT, &file_info,
                         err);
  if (file_id < 0 || !read_mdt_features(file_id, mdt, err)
      || !check_mdt_features(mdt, mlib, err)
      || !read_mdt_data(file_id, mdt, err)
      || !mdt_hdf_close(file_id, &file_info, err)) {
    return FALSE;
  } else {
    int ierr;
    mod_mdt_setup_check(mdt, &mlib->base, &ierr);
    if (ierr == 0) {
      return TRUE;
    } else {
      handle_modeller_error(err);
      return FALSE;
    }
  }
}
